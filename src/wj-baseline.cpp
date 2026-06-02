#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <zlib.h>
#include "kseq.h"

// 初始化 kseq，用于读取 FASTA/FASTQ 文件
KSEQ_INIT(gzFile, gzread)

// 并查集结构
class UnionFind {
public:
    std::vector<int> parent;

    UnionFind(int n) {
        parent.resize(n);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    // 查找并进行路径压缩
    int find(int i) {
        if (parent[i] == i) {
            return i;
        }
        return parent[i] = find(parent[i]); // 路径压缩
    }

    // 合并两个集合，强制根节点为较小的下标
    void unite(int i, int j) {
        int root_i = find(i);
        int root_j = find(j);
        if (root_i != root_j) {
            if (root_i < root_j) {
                parent[root_j] = root_i;
            } else {
                parent[root_i] = root_j;
            }
        }
    }

    // 最后执行一次全面的查找，将所有节点直接指向最终的根节点
    void flatten() {
        for (size_t i = 0; i < parent.size(); ++i) {
            find(i);
        }
    }
};

// 计算两个序列基于 k-mer 的 Weighted Jaccard 相似度
// 公式: sum(min(A_i, B_i)) / sum(max(A_i, B_i))
double weighted_jaccard(const std::string& s1, const std::string& s2, int k = 3) {
    if (s1.length() < k || s2.length() < k) return 0.0;

    std::unordered_map<std::string, int> kmers1;
    std::unordered_map<std::string, int> kmers2;

    for (size_t i = 0; i <= s1.length() - k; ++i) {
        kmers1[s1.substr(i, k)]++;
    }
    for (size_t i = 0; i <= s2.length() - k; ++i) {
        kmers2[s2.substr(i, k)]++;
    }

    long long min_sum = 0;
    long long max_sum = 0;

    for (const auto& pair : kmers1) {
        const std::string& kmer = pair.first;
        int count1 = pair.second;
        int count2 = kmers2[kmer]; // 如果不存在则为0
        min_sum += std::min(count1, count2);
        max_sum += std::max(count1, count2);
    }

    for (const auto& pair : kmers2) {
        if (kmers1.find(pair.first) == kmers1.end()) {
            max_sum += pair.second;
        }
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
    int kmer_size = 3; // 针对蛋白质通常用 3-mer

    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return 1;
    }

    kseq_t *seq = kseq_init(fp);
    std::vector<std::string> sequences;

    // 读取所有序列
    while (kseq_read(seq) >= 0) {
        sequences.push_back(seq->seq.s);
    }
    kseq_destroy(seq);
    gzclose(fp);

    int n = sequences.size();
    if (n == 0) {
        std::cerr << "No sequences found." << std::endl;
        return 0;
    }

    UnionFind uf(n);

    // O(N^2) 遍历所有序列对
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double sim = weighted_jaccard(sequences[i], sequences[j], kmer_size);
            if (sim >= threshold) {
                uf.unite(i, j);
		//debug only!
	    	//std::cerr << i << ":" << j << " " << sim << std::endl;
            }
        }
    }

    // 拍平并查集，确保所有路径被完全压缩到根
    uf.flatten();

    // 输出内部数组，用空格分隔
    for (int i = 0; i < n; ++i) {
        std::cout << uf.parent[i] << (i == n - 1 ? "" : " ");
    }
    std::cout << std::endl;

    return 0;
}
