/**
 * Weighted Jaccard Similarity Clustering - Phase 1 Deterministic
 * ==============================================================
 * This version keeps only deterministic optimizations:
 *   1. Direct k-mer state mapping for uppercase A-Z characters.
 *   2. 64-byte aligned contiguous frequency-vector storage.
 *   3. SIMD-friendly exact Weighted Jaccard accumulation.
 *   4. Cache-blocked O(N^2) pairwise comparison.
 *   5. Thread-local edge collection with deterministic union by smaller root.
 *
 * There is no probabilistic pre-filter in this version. Every sequence pair is
 * checked by exact Weighted Jaccard before union, so the result matches the
 * serial baseline for the same input and threshold.
 *
 * Compile: g++ -O3 -mavx2 -fopenmp -pthread wj.cpp -lz -o jaccard_cluster
 * Run:     ./jaccard_cluster <fasta_file> <threshold>
 * Or:      ./jaccard_cluster <number_of_threads> <fasta_file> <threshold>
 */

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

static constexpr int AA_NUM = 26;
static constexpr int KMER_DIM = AA_NUM * AA_NUM * AA_NUM;
static constexpr int K = 3;
static constexpr int BLOCK = 32;

static int aa_map[256];

static void init_aa_map() {
    std::memset(aa_map, -1, sizeof(aa_map));
    const char* letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    for (int i = 0; i < AA_NUM; ++i) {
        aa_map[static_cast<unsigned char>(letters[i])] = i;
    }
}

static inline int kmer_index(const char* s) {
    int c0 = aa_map[static_cast<unsigned char>(s[0])];
    int c1 = aa_map[static_cast<unsigned char>(s[1])];
    int c2 = aa_map[static_cast<unsigned char>(s[2])];
    if (c0 < 0 || c1 < 0 || c2 < 0) return -1;
    return c0 * AA_NUM * AA_NUM + c1 * AA_NUM + c2;
}

static void seq_to_freq(const std::string& seq, uint16_t* freq) {
    std::memset(freq, 0, KMER_DIM * sizeof(uint16_t));
    if (static_cast<int>(seq.size()) < K) return;

    const char* data = seq.data();
    const size_t limit = seq.size() - K;
    for (size_t i = 0; i <= limit; ++i) {
        int idx = kmer_index(data + i);
        if (idx >= 0) {
            ++freq[idx];
        }
    }
}

static inline double weighted_jaccard(const uint16_t* __restrict__ a,
                                      const uint16_t* __restrict__ b) {
    long long min_sum = 0;
    long long max_sum = 0;

    #pragma omp simd reduction(+:min_sum, max_sum)
    for (int i = 0; i < KMER_DIM; ++i) {
        uint16_t ai = a[i];
        uint16_t bi = b[i];
        min_sum += (ai < bi) ? ai : bi;
        max_sum += (ai > bi) ? ai : bi;
    }

    if (max_sum == 0) return 0.0;
    return static_cast<double>(min_sum) / static_cast<double>(max_sum);
}

class UnionFind {
public:
    explicit UnionFind(int n) : parent(n) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    int find(int i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }

    void unite(int i, int j) {
        int ri = find(i);
        int rj = find(j);
        if (ri == rj) return;
        if (ri < rj) {
            parent[rj] = ri;
        } else {
            parent[ri] = rj;
        }
    }

    void flatten() {
        for (size_t i = 0; i < parent.size(); ++i) {
            find(static_cast<int>(i));
        }
    }

    std::vector<int> parent;
};

static void print_usage(const char* program) {
    std::cerr << "Usage: " << program << " <fasta_file> <threshold>\n";
    std::cerr << "Alternative Usage: " << program
              << " <number_of_threads> <fasta_file> <threshold>\n";
}

int main(int argc, char* argv[]) {
    if (argc != 3 && argc != 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string filename;
    double threshold = 0.0;

    if (argc == 3) {
        filename = argv[1];
        threshold = std::stod(argv[2]);
    } else {
        int num_threads = std::stoi(argv[1]);
        if (num_threads <= 0) {
            std::cerr << "Thread count must be positive." << std::endl;
            return 1;
        }
        omp_set_num_threads(num_threads);
        filename = argv[2];
        threshold = std::stod(argv[3]);
    }

    init_aa_map();

    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return 1;
    }

    kseq_t* seq = kseq_init(fp);
    std::vector<std::string> sequences;
    while (kseq_read(seq) >= 0) {
        sequences.emplace_back(seq->seq.s);
    }
    kseq_destroy(seq);
    gzclose(fp);

    const int n = static_cast<int>(sequences.size());
    if (n == 0) {
        std::cerr << "No sequences found." << std::endl;
        return 0;
    }

    const size_t row_bytes = static_cast<size_t>(KMER_DIM) * sizeof(uint16_t);
    const size_t padded_row = (row_bytes + 63) & ~static_cast<size_t>(63);
    const size_t padded_elems = padded_row / sizeof(uint16_t);
    const size_t total_bytes = static_cast<size_t>(n) * padded_row;

    uint16_t* freq_pool = nullptr;
#ifdef _WIN32
    freq_pool = static_cast<uint16_t*>(_aligned_malloc(total_bytes, 64));
    if (!freq_pool) {
        std::cerr << "Memory allocation failed." << std::endl;
        return 1;
    }
#else
    {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, 64, total_bytes) != 0 || !ptr) {
            std::cerr << "Memory allocation failed." << std::endl;
            return 1;
        }
        freq_pool = static_cast<uint16_t*>(ptr);
    }
#endif

    std::memset(freq_pool, 0, total_bytes);

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i) {
        seq_to_freq(sequences[i],
                    freq_pool + static_cast<size_t>(i) * padded_elems);
    }

    sequences.clear();
    sequences.shrink_to_fit();

    UnionFind uf(n);

    #pragma omp parallel
    {
        std::vector<std::pair<int, int>> local_edges;

        #pragma omp for schedule(dynamic) nowait
        for (int bi = 0; bi < n; bi += BLOCK) {
            for (int bj = bi; bj < n; bj += BLOCK) {
                const int i_end = std::min(bi + BLOCK, n);
                const int j_end = std::min(bj + BLOCK, n);
                for (int i = bi; i < i_end; ++i) {
                    const uint16_t* vec_i =
                        freq_pool + static_cast<size_t>(i) * padded_elems;
                    const int j_start = (bi == bj) ? i + 1 : bj;
                    for (int j = j_start; j < j_end; ++j) {
                        const uint16_t* vec_j =
                            freq_pool + static_cast<size_t>(j) * padded_elems;
                        double sim = weighted_jaccard(vec_i, vec_j);
                        if (sim >= threshold) {
                            local_edges.emplace_back(i, j);
                        }
                    }
                }
            }
        }

        #pragma omp critical
        {
            for (const auto& edge : local_edges) {
                uf.unite(edge.first, edge.second);
            }
        }
    }

    uf.flatten();

    for (int i = 0; i < n; ++i) {
        std::cout << uf.parent[i] << (i == n - 1 ? "" : " ");
    }
    std::cout << std::endl;

#ifdef _WIN32
    _aligned_free(freq_pool);
#else
    std::free(freq_pool);
#endif

    return 0;
}
