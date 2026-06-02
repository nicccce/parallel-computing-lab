#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
// #include "./zlib-1.3.2/zlib.h"
#include <zlib.h>
#include "kseq.h"
#include <omp.h>

KSEQ_INIT(gzFile, gzread)

// ---------- 氨基酸映射 ----------
constexpr int AA_NUM = 20;
constexpr int KMER_DIM = AA_NUM * AA_NUM * AA_NUM;  // 8000

int aa_map[256];
void init_aa_map() {
    memset(aa_map, -1, sizeof(aa_map));
    const char* letters = "ACDEFGHIKLMNPQRSTVWY";
    for (int i = 0; i < AA_NUM; ++i)
        aa_map[(unsigned char)letters[i]] = i;
}

inline int kmer_index(const char* s) {
    return aa_map[(unsigned char)s[0]] * 400 +
           aa_map[(unsigned char)s[1]] * 20 +
           aa_map[(unsigned char)s[2]];
}

// 将序列转为频率向量（长度必须 >= k，不满足的序列向量全 0）
void seq_to_freq(const std::string& seq, uint16_t* freq, int k = 3) {
    memset(freq, 0, KMER_DIM * sizeof(uint16_t));
    if (seq.size() < k) return;
    for (size_t i = 0; i <= seq.size() - k; ++i) {
        int idx = kmer_index(&seq[i]);
        if (idx >= 0) freq[idx]++;   // 忽略含非法字符的 k-mer
    }
}

// ---------- 并查集 ----------
class UnionFind {
public:
    std::vector<int> parent;
    UnionFind(int n) : parent(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }
    void unite(int i, int j) {
        int ri = find(i), rj = find(j);
        if (ri != rj) {
            if (ri < rj) parent[rj] = ri;
            else         parent[ri] = rj;
        }
    }
    void flatten() {
        for (size_t i = 0; i < parent.size(); ++i)
            find(i);
    }
};

// ---------- 向量化相似度 ----------
inline double weighted_jaccard(const uint16_t* a, const uint16_t* b) {
    long long min_sum = 0, max_sum = 0;
    #pragma omp simd reduction(+:min_sum, max_sum)
    for (int i = 0; i < KMER_DIM; ++i) {
        uint16_t ai = a[i], bi = b[i];
        min_sum += ai < bi ? ai : bi;
        max_sum += ai > bi ? ai : bi;
    }
    if (max_sum == 0) return 0.0;
    return static_cast<double>(min_sum) / max_sum;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <threshold>" << std::endl;
        return 1;
    }
    std::string filename = argv[1];
    double threshold = std::stod(argv[2]);

    init_aa_map();

    // ---------- 读取序列 ----------
    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) { std::cerr << "Failed to open file\n"; return 1; }
    kseq_t *seq = kseq_init(fp);
    std::vector<std::string> sequences;
    while (kseq_read(seq) >= 0)
        sequences.emplace_back(seq->seq.s);
    kseq_destroy(seq);
    gzclose(fp);

    int n = sequences.size();
    if (n == 0) { std::cerr << "No sequences.\n"; return 0; }

    // ---------- 预计算频率矩阵（行主序）----------
    std::vector<uint16_t> freqs(n * KMER_DIM, 0);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i)
        seq_to_freq(sequences[i], &freqs[i * KMER_DIM]);

    // 释放原始序列内存（可选）
    sequences.clear();
    sequences.shrink_to_fit();

    // ---------- 缓存分块并行比较 ----------
    constexpr int BLOCK_SIZE = 512;   // 根据 L3 缓存调整
    UnionFind uf(n);

    // 边列表（每个线程局部收集，避免竞争）
    std::vector<std::pair<int,int>> all_edges;
    #pragma omp parallel
    {
        std::vector<std::pair<int,int>> local_edges;
        #pragma omp for schedule(dynamic) nowait
        for (int bi = 0; bi < n; bi += BLOCK_SIZE) {
            for (int bj = bi; bj < n; bj += BLOCK_SIZE) {
                int i_end = std::min(bi + BLOCK_SIZE, n);
                int j_end = std::min(bj + BLOCK_SIZE, n);
                for (int i = bi; i < i_end; ++i) {
                    const uint16_t* vec_i = &freqs[i * KMER_DIM];
                    int j_start = (bi == bj) ? i + 1 : bj;
                    for (int j = j_start; j < j_end; ++j) {
                        double sim = weighted_jaccard(vec_i, &freqs[j * KMER_DIM]);
                        if (sim >= threshold)
                            local_edges.emplace_back(i, j);
                    }
                }
            }
        }
        // 合并到全局边列表（临界区）
        #pragma omp critical
        all_edges.insert(all_edges.end(), local_edges.begin(), local_edges.end());
    }

    // 串行合并并查集
    for (auto& e : all_edges)
        uf.unite(e.first, e.second);

    uf.flatten();

    // 输出
    for (int i = 0; i < n; ++i)
        std::cout << uf.parent[i] << (i == n-1 ? "" : " ");
    std::cout << std::endl;

    return 0;
}