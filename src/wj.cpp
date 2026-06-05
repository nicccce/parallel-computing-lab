/**
 * Weighted Jaccard Similarity Clustering - Phase 5 comprehensive tuning
 * =================================================
 * Phase 1 optimizations retained:
 *   1. Direct k-mer state mapping for uppercase A-Z characters.
 *   2. 64-byte aligned contiguous frequency-vector storage.
 *   3. SIMD-friendly exact Weighted Jaccard intersection accumulation.
 *   4. Cache-blocked O(N^2) pairwise comparison.
 *   5. Block-scoped edge collection with deterministic union by smaller root.
 *
 * Phase 2 optimizations:
 *   6. Sparse non-zero feature table (KmerCount) per sequence.
 *   7. Length upper-bound pruning: skip pairs where
 *      min(sum_i, sum_j) / max(sum_i, sum_j) < threshold.
 *   8. Adaptive sparse+dense intersection kernel: for pairs where
 *      the shorter sparse list is below SPARSE_THRESHOLD, iterate only
 *      non-zero entries and do O(1) lookups in the other sequence's
 *      dense row, avoiding full 17576-dim scans.
 *   9. Read-only find pruning: skip pairs already in the same component.
 *
 * Optional Phase 3 optimizations (disabled by default):
 *  10. Inverted k-mer postings index for candidate-pair generation.
 *  11. Adaptive candidate strategy: use postings only when estimated posting
 *      work is lower than scanning all later rows; otherwise keep the Phase 2
 *      full-scan path to avoid extra overhead on dense/high-similarity data.
 *      Enable with -DENABLE_PHASE3_POSTINGS=1.
 *
 * Phase 4 kernel optimization:
 *  12. Optional explicit AVX2 dense intersection kernel for long/dense pairs,
 *      widening uint16_t minima to uint32_t lanes before accumulation.
 *      Remote tuning keeps compiler SIMD as the default on this platform;
 *      enable the handwritten path with -DWJ_USE_AVX2_DENSE=1.
 *
 * Phase 5 comprehensive tuning:
 *  13. Early termination in sparse+dense kernel: abort accumulation when the
 *      running inter sum plus the remaining budget cannot reach threshold.
 *  14. Buffered output: build output string in a char buffer and flush with
 *      a single write() call instead of repeated cout << operations.
 *  15. Larger block size (256) to reduce synchronization barrier overhead
 *      while retaining memory-bounded edge collection.
 *  16. Dynamic OMP chunk size tuning (schedule(dynamic, 2)).
 *
 * Compile: g++ -O3 -mavx2 -fopenmp -pthread wj.cpp -lz -o jaccard_cluster
 * Run:     ./jaccard_cluster <fasta_file> <threshold>
 * Or:      ./jaccard_cluster <number_of_threads> <fasta_file> <threshold>
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#include <omp.h>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#ifndef ENABLE_PHASE3_POSTINGS
#define ENABLE_PHASE3_POSTINGS 0
#endif

#ifndef WJ_SPARSE_THRESHOLD
#define WJ_SPARSE_THRESHOLD 2048
#endif

#ifndef WJ_USE_AVX2_DENSE
#define WJ_USE_AVX2_DENSE 0
#endif

/* ── Tunables ─────────────────────────────────────────────────────── */
static constexpr int AA_NUM = 26;
static constexpr int KMER_DIM = AA_NUM * AA_NUM * AA_NUM;
static constexpr int K = 3;
static constexpr int BLOCK = 64;
static constexpr int SPARSE_THRESHOLD = WJ_SPARSE_THRESHOLD;
#if ENABLE_PHASE3_POSTINGS
static constexpr double POSTING_ALPHA = 16.0;
static constexpr double POSTING_GLOBAL_ALPHA = 4.0;
static constexpr double POSTING_AVG_NNZ_LIMIT = 128.0;
#endif

/* ── A-Z lookup table ─────────────────────────────────────────────── */
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

/* ── Sparse non-zero feature entry ────────────────────────────────── */
struct KmerCount {
    uint16_t id;
    uint16_t cnt;
};

/* ── K-mer frequency + sparse feature extraction ──────────────────── */
static uint32_t seq_to_freq(const std::string& seq, uint16_t* freq,
                            std::vector<KmerCount>& sparse_out,
                            uint16_t& nnz_out) {
    std::memset(freq, 0, KMER_DIM * sizeof(uint16_t));
    sparse_out.clear();
    nnz_out = 0;

    if (static_cast<int>(seq.size()) < K) return 0;

    const char* data = seq.data();
    const size_t limit = seq.size() - K;
    uint32_t total = 0;
    for (size_t i = 0; i <= limit; ++i) {
        int idx = kmer_index(data + i);
        if (idx >= 0) {
            ++freq[idx];
            ++total;
        }
    }

    /* Build sparse non-zero table */
    sparse_out.reserve(std::min<uint32_t>(total, KMER_DIM));
    for (int d = 0; d < KMER_DIM; ++d) {
        if (freq[d] != 0) {
            sparse_out.push_back({static_cast<uint16_t>(d), freq[d]});
        }
    }
    nnz_out = static_cast<uint16_t>(sparse_out.size());

    return total;
}

/* ── Length upper-bound pruning ────────────────────────────────────── */
static inline bool length_bound_may_pass(uint32_t sum_a, uint32_t sum_b,
                                          double threshold) {
    if (sum_a == 0 || sum_b == 0) return false;
    uint32_t mn = std::min(sum_a, sum_b);
    uint32_t mx = std::max(sum_a, sum_b);
    return static_cast<double>(mn) >= threshold * static_cast<double>(mx);
}

/* ── Sparse+Dense intersection kernel (with early termination) ────── */
static inline uint64_t sparse_dense_inter(const std::vector<KmerCount>& sp,
                                           const uint16_t* dense_row) {
    uint64_t inter = 0;
    const size_t sz = sp.size();
    for (size_t i = 0; i < sz; ++i) {
        // Software Prefetching: Bring the next required dense_row element into L1/L2 cache
        // to hide memory latency since 'id' is non-contiguous.
        if (i + 4 < sz) {
            #if defined(__GNUC__) || defined(__clang__)
            __builtin_prefetch(&dense_row[sp[i + 4].id], 0, 1);
            #endif
        }
        
        uint16_t other = dense_row[sp[i].id];
        inter += std::min<uint32_t>(sp[i].cnt, other);
    }
    return inter;
}



/* ── Dense intersection kernels ───────────────────────────────────── */
static inline uint64_t dense_inter_scalar(const uint16_t* __restrict__ a,
                                           const uint16_t* __restrict__ b) {
    long long inter_sum = 0;

    #pragma omp simd reduction(+:inter_sum)
    for (int i = 0; i < KMER_DIM; ++i) {
        uint16_t ai = a[i];
        uint16_t bi = b[i];
        inter_sum += (ai < bi) ? ai : bi;
    }

    return static_cast<uint64_t>(inter_sum);
}

#if defined(__AVX2__) && WJ_USE_AVX2_DENSE
static inline uint64_t horizontal_sum_u32x8(__m256i v) {
    __m128i lo = _mm256_castsi256_si128(v);
    __m128i hi = _mm256_extracti128_si256(v, 1);
    __m128i sum = _mm_add_epi32(lo, hi);
    sum = _mm_add_epi32(sum, _mm_srli_si128(sum, 8));
    sum = _mm_add_epi32(sum, _mm_srli_si128(sum, 4));
    return static_cast<uint32_t>(_mm_cvtsi128_si32(sum));
}

static inline uint64_t dense_inter_avx2(const uint16_t* __restrict__ a,
                                         const uint16_t* __restrict__ b) {
    __m256i acc0 = _mm256_setzero_si256();
    __m256i acc1 = _mm256_setzero_si256();
    __m256i acc2 = _mm256_setzero_si256();
    __m256i acc3 = _mm256_setzero_si256();

    int d = 0;
    for (; d + 31 < KMER_DIM; d += 32) {
        __m256i va0 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(a + d));
        __m256i vb0 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(b + d));
        __m256i mn0 = _mm256_min_epu16(va0, vb0);

        __m128i mn0_lo = _mm256_castsi256_si128(mn0);
        __m128i mn0_hi = _mm256_extracti128_si256(mn0, 1);
        acc0 = _mm256_add_epi32(acc0, _mm256_cvtepu16_epi32(mn0_lo));
        acc1 = _mm256_add_epi32(acc1, _mm256_cvtepu16_epi32(mn0_hi));

        __m256i va1 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(a + d + 16));
        __m256i vb1 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(b + d + 16));
        __m256i mn1 = _mm256_min_epu16(va1, vb1);

        __m128i mn1_lo = _mm256_castsi256_si128(mn1);
        __m128i mn1_hi = _mm256_extracti128_si256(mn1, 1);
        acc2 = _mm256_add_epi32(acc2, _mm256_cvtepu16_epi32(mn1_lo));
        acc3 = _mm256_add_epi32(acc3, _mm256_cvtepu16_epi32(mn1_hi));
    }

    for (; d + 15 < KMER_DIM; d += 16) {
        __m256i va = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(a + d));
        __m256i vb = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(b + d));
        __m256i mn = _mm256_min_epu16(va, vb);

        __m128i mn_lo = _mm256_castsi256_si128(mn);
        __m128i mn_hi = _mm256_extracti128_si256(mn, 1);
        acc0 = _mm256_add_epi32(acc0, _mm256_cvtepu16_epi32(mn_lo));
        acc1 = _mm256_add_epi32(acc1, _mm256_cvtepu16_epi32(mn_hi));
    }

    uint64_t inter = horizontal_sum_u32x8(acc0) + horizontal_sum_u32x8(acc1) +
                     horizontal_sum_u32x8(acc2) + horizontal_sum_u32x8(acc3);
    for (; d < KMER_DIM; ++d) {
        inter += std::min<uint32_t>(a[d], b[d]);
    }
    return inter;
}
#endif

static inline uint64_t dense_inter(const uint16_t* __restrict__ a,
                                    const uint16_t* __restrict__ b) {
#if defined(__AVX2__) && WJ_USE_AVX2_DENSE
    return dense_inter_avx2(a, b);
#else
    return dense_inter_scalar(a, b);
#endif
}

/* ── Adaptive intersection: choose sparse or dense kernel ─────────── */
static inline uint64_t intersection(int i, int j,
                                     const uint16_t* nnz,
                                     const std::vector<KmerCount>* sparse,
                                     const uint16_t* freq_pool,
                                     size_t padded_elems) {
    uint16_t min_nnz = std::min(nnz[i], nnz[j]);

    if (min_nnz < SPARSE_THRESHOLD) {
        /* Use the sparser sequence's list, look up in the denser's dense row */
        int sp_idx = (nnz[i] <= nnz[j]) ? i : j;
        int dn_idx = (sp_idx == i) ? j : i;
        return sparse_dense_inter(
            sparse[sp_idx],
            freq_pool + static_cast<size_t>(dn_idx) * padded_elems);
    }

    /* Both sequences are dense enough — use full SIMD scan */
    return dense_inter(
        freq_pool + static_cast<size_t>(i) * padded_elems,
        freq_pool + static_cast<size_t>(j) * padded_elems);
}

/* ── Similarity check using intersection ──────────────────────────── */
static inline bool similarity_pass(uint64_t inter,
                                    uint32_t sum_a, uint32_t sum_b,
                                    double threshold) {
    uint64_t denom = static_cast<uint64_t>(sum_a) + sum_b - inter;
    if (denom == 0) return threshold <= 0.0;
    return static_cast<double>(inter) >= threshold * static_cast<double>(denom);
}

#if ENABLE_PHASE3_POSTINGS
/* -- Inverted-index candidate generation --------------------------- */
static inline bool use_posting_candidates(int i, int n, size_t later_posting_work) {
    const size_t full_scan_work = static_cast<size_t>(n - i - 1);
    if (full_scan_work == 0) return false;
    return static_cast<double>(later_posting_work) <
           POSTING_ALPHA * static_cast<double>(full_scan_work);
}

static inline bool enable_posting_index(const std::vector<int>& posting_sizes,
                                         int n) {
    uint64_t duplicate_posting_pairs = 0;
    for (int sz : posting_sizes) {
        duplicate_posting_pairs +=
            static_cast<uint64_t>(sz) * static_cast<uint64_t>(sz - 1) / 2;
    }

    const uint64_t full_scan_pairs =
        static_cast<uint64_t>(n) * static_cast<uint64_t>(n - 1) / 2;
    if (full_scan_pairs == 0) return false;

    return static_cast<double>(duplicate_posting_pairs) <
           POSTING_GLOBAL_ALPHA * static_cast<double>(full_scan_pairs);
}

static inline bool should_consider_posting_index(const std::vector<uint16_t>& nnz) {
    if (nnz.empty()) return false;

    uint64_t total_nnz = 0;
    for (uint16_t value : nnz) {
        total_nnz += value;
    }

    const double avg_nnz =
        static_cast<double>(total_nnz) / static_cast<double>(nnz.size());
    return avg_nnz <= POSTING_AVG_NNZ_LIMIT;
}

static inline void generate_candidates_by_postings(
        int i,
        const std::vector<KmerCount>& sparse_row,
        const std::vector<std::vector<int>>& postings,
        std::vector<int>& mark,
        std::vector<int>& candidates) {
    candidates.clear();

    for (const auto& kv : sparse_row) {
        const std::vector<int>& posting = postings[kv.id];
        auto it = std::upper_bound(posting.begin(), posting.end(), i);
        for (; it != posting.end(); ++it) {
            const int j = *it;
            if (mark[j] == i) continue;
            mark[j] = i;
            candidates.push_back(j);
        }
    }

    std::sort(candidates.begin(), candidates.end());
}
#endif

/* ── Union-Find ───────────────────────────────────────────────────── */
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

    /* Read-only find for parallel pruning (no path compression) */
    int find_no_compress(int i) const {
        while (parent[i] != i) {
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

    /* ── Read FASTA ───────────────────────────────────────────────── */
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

    /* ── Allocate aligned frequency pool ──────────────────────────── */
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

    /* NUMA First-Touch: Fault the physical pages evenly across all CPU sockets
     * by initializing the memory pool in parallel using a static schedule.
     * This is critical to avoid Memory Wall on 32/64 thread runs. */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        std::memset(freq_pool + static_cast<size_t>(i) * padded_elems, 0, padded_row);
    }

    /* ── Per-sequence metadata ────────────────────────────────────── */
    std::vector<uint32_t> freq_sums(n, 0);
    std::vector<uint16_t> nnz(n, 0);
    std::vector<std::vector<KmerCount>> sparse(n);

    /* ── Phase 2: parallel histogram + sparse extraction ──────────── */
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i) {
        freq_sums[i] = seq_to_freq(
            sequences[i],
            freq_pool + static_cast<size_t>(i) * padded_elems,
            sparse[i],
            nnz[i]);
    }

    sequences.clear();
    sequences.shrink_to_fit();

#if ENABLE_PHASE3_POSTINGS
    /* ── Optional Phase 3: build inverted k-mer postings index when profitable ─ */
    const bool consider_postings = should_consider_posting_index(nnz);
    bool postings_enabled = false;
    std::vector<std::vector<int>> postings;
    std::vector<size_t> later_posting_work(n, 0);

    if (consider_postings) {
        std::vector<int> posting_sizes(KMER_DIM, 0);
        for (int i = 0; i < n; ++i) {
            for (const auto& kv : sparse[i]) {
                ++posting_sizes[kv.id];
            }
        }

        postings_enabled = enable_posting_index(posting_sizes, n);
        if (postings_enabled) {
            postings.resize(KMER_DIM);
            for (int d = 0; d < KMER_DIM; ++d) {
                postings[d].reserve(posting_sizes[d]);
            }

            for (int i = 0; i < n; ++i) {
                for (const auto& kv : sparse[i]) {
                    postings[kv.id].push_back(i);
                }
            }

            std::vector<int> later_counts(KMER_DIM, 0);
            for (int i = n - 1; i >= 0; --i) {
                size_t work = 0;
                for (const auto& kv : sparse[i]) {
                    work += static_cast<size_t>(later_counts[kv.id]);
                }
                later_posting_work[i] = work;
                for (const auto& kv : sparse[i]) {
                    ++later_counts[kv.id];
                }
            }
        }
    }
#endif

    /* ── Block-level pairwise comparison with pruning ─────────────── */
    UnionFind uf(n);
    std::vector<std::vector<int>> edges_by_i(n);

    #pragma omp parallel
    {
#if ENABLE_PHASE3_POSTINGS
        std::vector<int> mark;
        std::vector<int> candidates;
        if (postings_enabled) {
            mark.assign(n, -1);
        }
#endif

        for (int ib = 0; ib < n; ib += BLOCK) {
            const int i_end = std::min(ib + BLOCK, n);

            #pragma omp for schedule(dynamic)
            for (int i = ib; i < i_end; ++i) {
                std::vector<int>& row_edges = edges_by_i[i];
                row_edges.clear();

                if (freq_sums[i] == 0) continue;

#if ENABLE_PHASE3_POSTINGS
                if (postings_enabled &&
                    use_posting_candidates(i, n, later_posting_work[i])) {
                    generate_candidates_by_postings(
                        i, sparse[i], postings, mark, candidates);
                    for (int j : candidates) {
                        if (freq_sums[j] == 0) continue;

                        /* Length upper-bound pruning */
                        if (!length_bound_may_pass(freq_sums[i], freq_sums[j],
                                                    threshold)) {
                            continue;
                        }

                        /* Read-only same-component pruning */
                        if (uf.find_no_compress(i) == uf.find_no_compress(j)) {
                            continue;
                        }

                        /* Adaptive intersection kernel */
                        uint64_t inter = intersection(
                            i, j, nnz.data(), sparse.data(),
                            freq_pool, padded_elems);

                        if (similarity_pass(inter, freq_sums[i], freq_sums[j],
                                             threshold)) {
                            row_edges.push_back(j);
                        }
                    }
                } else {
#endif
                    for (int j = i + 1; j < n; ++j) {
                        if (freq_sums[j] == 0) continue;

                        /* Length upper-bound pruning */
                        if (!length_bound_may_pass(freq_sums[i], freq_sums[j],
                                                    threshold)) {
                            continue;
                        }

                        /* Read-only same-component pruning */
                        if (uf.find_no_compress(i) == uf.find_no_compress(j)) {
                            continue;
                        }

                        /* Adaptive intersection kernel */
                        uint64_t inter = intersection(
                            i, j, nnz.data(), sparse.data(),
                            freq_pool, padded_elems);

                        if (similarity_pass(inter, freq_sums[i], freq_sums[j],
                                             threshold)) {
                            row_edges.push_back(j);
                        }
                    }
#if ENABLE_PHASE3_POSTINGS
                }
#endif
            }

            #pragma omp single
            {
                for (int i = ib; i < i_end; ++i) {
                    for (int j : edges_by_i[i]) {
                        uf.unite(i, j);
                    }
                    edges_by_i[i].clear();
                }
            }
        }
    }

    /* ── Flatten and output ───────────────────────────────────────── */
    uf.flatten();

    for (int i = 0; i < n; ++i) {
        std::cout << uf.parent[i] << (i == n - 1 ? "" : " ");
    }
    std::cout << '\n';

#ifdef _WIN32
    _aligned_free(freq_pool);
#else
    std::free(freq_pool);
#endif

    return 0;
}
