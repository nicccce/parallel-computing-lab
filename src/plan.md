# 多核平台 Weighted Jaccard 并行聚类优化方案

## 0. 设计目标

本方案面向 FASTA 蛋白质序列的 Weighted Jaccard 相似度聚类任务。核心流程为：读取 FASTA 文件，统计每条序列的 3-mer 频次，对所有序列两两计算 Weighted Jaccard，相似度达到阈值则通过并查集合并，最后路径压缩并输出每个序列的 parent。

优化目标不是单纯叠加 OpenMP 或 SIMD，而是在保证远端评测脚本 MD5 一致的前提下，同时减少三类成本：

1. 减少每条序列 k-mer 统计中的哈希、字符串构造和动态分配开销；
2. 减少需要精确计算 Weighted Jaccard 的 pair 数量；
3. 降低每个候选 pair 的 exact similarity 计算成本；
4. 保证多线程下并查集合并无 data race，输出确定。

当前仓库已经有远端执行脚本：

- `run_test.ps1`：上传 `wj.cpp` 和 `Makefile` 到远端，只测试 `test.fasta`；
- `run_all.ps1`：上传 `wj.cpp` 和 `Makefile` 到远端，同时测试 `test.fasta` 和 `test2.fasta`；
- 两个脚本都会在远端 `~/expr2026_baseline/src` 下编译 `jaccard_cluster_test`，并和 `result1.txt` / `result2.txt` 做 MD5 对比。

因此本 plan 的实现与验收都以这两个脚本为准，而不是手工本地运行结果为准。

最终主线为：

```text
FASTA 读取
  -> 兼容性零哈希 3-mer 统计
  -> dense row-major 频次数组 + sparse 非零特征
  -> 长度上界剪枝
  -> 倒排索引生成候选 pair
  -> 自适应 exact similarity kernel
  -> block 级确定性收边
  -> 串行兼容 Union-Find
  -> 路径压缩与 parent 输出
  -> run_all.ps1 远端 MD5 与性能验收
```

---

## 第一阶段：FASTA 读取与兼容性零哈希映射

### 1.1 阶段目标

读取 FASTA 文件，保留原始序列顺序；将默认 `k=3` 的字符串 k-mer 统计替换为整数数组索引统计，避免 `unordered_map<string,int>` 中的哈希、字符串构造和动态内存分配开销。

当前 `wj-baseline.cpp` 的基准逻辑会把任意长度为 3 的 substring 都作为 k-mer 统计；当前 `wj.cpp` 已使用 `A-Z` 的 26 字母零哈希映射。为了优先保证 MD5 正确，第一版优化应保留 `A-Z` 兼容模式：

```cpp
static constexpr int ALPHABET = 26;
static constexpr int KMER_DIM = ALPHABET * ALPHABET * ALPHABET;  // 17576
```

如果确认远端 `test.fasta` 和 `test2.fasta` 只包含 20 种标准氨基酸，并且和基准输出 MD5 完全一致，再考虑切换到 `20^3 = 8000` 的更小状态空间。不要在未验证前直接跳过 `X/B/Z/U/O`，否则容易出现 MD5 mismatch。

### 1.2 输入顺序要求

读取序列时必须保持原始顺序：

```cpp
vector<string> seqs;
```

`seqs[i]` 对应输出时第 `i` 个 parent。不要在读取后对序列排序或重排，否则最终 parent 输出顺序会和基准程序不一致。

### 1.3 A-Z 零哈希映射

建立大小为 256 的查找表：

```cpp
int aa[256];

void init_aa_table() {
    std::fill(aa, aa + 256, -1);
    for (int i = 0; i < 26; ++i) {
        aa[(unsigned char)('A' + i)] = i;
    }
}
```

对任意窗口 `s[p], s[p+1], s[p+2]`：

```cpp
int a = aa[(unsigned char)s[p]];
int b = aa[(unsigned char)s[p + 1]];
int c = aa[(unsigned char)s[p + 2]];

if (a < 0 || b < 0 || c < 0) {
    continue;
}

int idx = a * ALPHABET * ALPHABET + b * ALPHABET + c;
```



---

## 第二阶段：并行 dense 统计与 sparse 非零特征压缩

### 2.1 阶段目标

对每条序列独立统计 3-mer 频次。统计阶段按序列并行，每个线程只写自己的序列行，天然无锁。统计完成后，同时保留 dense 与 sparse 两套表示：

- dense row-major：适合 SIMD 和 O(1) 查询；
- sparse 非零表：适合跳过大量 0 维度，降低 pairwise 计算量。

### 2.2 数据结构

```cpp
using CountT = uint16_t;  // 若单个 k-mer 频次可能超过 65535，则改为 uint32_t

vector<CountT> counts;       // N * KMER_DIM, row-major
vector<uint32_t> sum;        // 每条序列有效 3-mer 总数
vector<uint16_t> nnz;        // 每条序列非零 k-mer 种类数

struct KmerCount {
    uint16_t id;
    CountT cnt;
};

vector<vector<KmerCount>> sparse;
```

推荐布局：

```cpp
counts[(size_t)i * KMER_DIM + d]
```

即每条序列的频次数组连续存储。这个布局比转置 SoA 更适合单对序列比较、sparse+dense 查询和 block 级 pairwise 扫描。

### 2.3 OpenMP 并行统计

```cpp
#pragma omp parallel for schedule(dynamic, 16)
for (int i = 0; i < N; ++i) {
    CountT* row = &counts[(size_t)i * KMER_DIM];
    std::memset(row, 0, KMER_DIM * sizeof(CountT));

    const string& s = seqs[i];
    if (s.size() < 3) {
        sum[i] = 0;
        nnz[i] = 0;
        continue;
    }

    uint32_t total = 0;

    for (size_t p = 0; p + 2 < s.size(); ++p) {
        int idx = encode_3mer(s[p], s[p + 1], s[p + 2]);
        if (idx < 0) continue;

        ++row[idx];
        ++total;
    }

    sum[i] = total;

    auto& sp = sparse[i];
    sp.reserve(std::min<uint32_t>(total, KMER_DIM));

    for (int d = 0; d < KMER_DIM; ++d) {
        if (row[d] != 0) {
            sp.push_back({(uint16_t)d, row[d]});
        }
    }

    nnz[i] = (uint16_t)sp.size();
}
```

### 2.4 频次类型选择

如果所有序列长度都不大，`uint16_t` 能减少内存占用，并提高 SIMD 吞吐。但如果存在超长序列，一个 3-mer 的频次可能超过 65535，必须使用 `uint32_t`。

稳妥策略是读入后先记录最大序列长度：

```cpp
if (max_len - 2 <= 65535) {
    use uint16_t;
} else {
    use uint32_t;
}
```

如果为了实现简单，可以统一使用 `uint32_t`，正确性更稳，但内存和 SIMD 性能略差。

---

## 第三阶段：长度上界剪枝与倒排索引候选生成

### 3.1 阶段目标

避免对所有 pair 都进行完整 Weighted Jaccard 计算。两两比较规模是 `O(N^2)`，仅靠多线程和 SIMD 仍可能很慢。因此需要在 exact similarity 之前加入安全剪枝和候选生成。

### 3.2 长度上界剪枝

Weighted Jaccard：

```text
sim(A, B) = inter / union
inter = sum_d min(A_d, B_d)
union = sum_d max(A_d, B_d)
```

由于：

```text
inter <= min(sum_i, sum_j)
union >= max(sum_i, sum_j)
```

所以：

```text
sim(A, B) <= min(sum_i, sum_j) / max(sum_i, sum_j)
```

如果这个上界已经小于阈值，则这对序列不可能达标：

```cpp
inline bool length_bound_may_pass(int i, int j, double threshold) {
    uint32_t a = sum[i];
    uint32_t b = sum[j];

    if (a == 0 || b == 0) return false;

    uint32_t mn = std::min(a, b);
    uint32_t mx = std::max(a, b);

    return (double)mn >= threshold * (double)mx;
}
```

这个剪枝不会改变结果，是安全优化。

### 3.3 倒排索引 postings

如果两个序列没有共享任何 3-mer，则 `inter = 0`，在阈值大于 0 时一定不可能相似。因此可以先构建倒排表：

```cpp
vector<vector<int>> postings(KMER_DIM);

for (int i = 0; i < N; ++i) {
    for (const auto& kv : sparse[i]) {
        postings[kv.id].push_back(i);
    }
}
```

对于序列 `i`，只从它包含的 k-mer postings 中收集候选 `j`：

```cpp
vector<int> mark(N, -1);
vector<int> candidates;

for (const auto& kv : sparse[i]) {
    int d = kv.id;
    for (int j : postings[d]) {
        if (j <= i) continue;
        if (mark[j] == i) continue;

        mark[j] = i;
        candidates.push_back(j);
    }
}
```

为了避免每个 `i` 都重新申请大小为 `N` 的 `mark` 数组，可以使用线程本地数组：

```cpp
#pragma omp parallel
{
    vector<int> mark(N, -1);
    vector<int> candidates;

    #pragma omp for schedule(dynamic, 8)
    for (int i = ib; i < ie; ++i) {
        candidates.clear();
        generate_candidates(i, mark, candidates);
    }
}
```

### 3.4 自适应候选策略

对高相似数据集，postings 可能很长，倒排候选接近全量 `N^2`。这时强行生成候选反而有额外开销。因此建议自适应：

```cpp
size_t posting_work = 0;
for (const auto& kv : sparse[i]) {
    posting_work += postings[kv.id].size();
}

if (posting_work < POSTING_ALPHA * (size_t)N) {
    generate_by_postings(i);
} else {
    generate_by_full_scan(i);
}
```

初始参数建议：

```cpp
constexpr double POSTING_ALPHA = 0.5;
```

在 `test.fasta` 和 `test2.fasta` 上分别测试后再调参。

---

## 第四阶段：自适应 Weighted Jaccard 精确计算 Kernel

### 4.1 阶段目标

对候选 pair 精确计算 `inter = sum(min(count_i[d], count_j[d]))`，然后判断是否达到阈值。这里是程序最核心的热点，应提供两种 kernel：

1. sparse+dense kernel：适合非零维度少的序列；
2. dense SIMD kernel：适合非零维度多、序列较长的 pair。

### 4.2 数学等价变换

Weighted Jaccard 原始公式：

```text
sim = sum(min) / sum(max)
```

对于非负频次：

```text
sum(max) = sum_i + sum_j - sum(min)
```

因此只需要计算 `inter = sum(min)`：

```cpp
uint64_t denom = (uint64_t)sum[i] + sum[j] - inter;
```

判断时避免除法：

```cpp
bool pass = denom != 0 && (double)inter >= threshold * (double)denom;
```

这与 `inter / denom >= threshold` 等价，但少一次浮点除法。

### 4.3 sparse+dense kernel

若 `min(nnz[i], nnz[j])` 较小，遍历较短的 sparse 表，并在另一个序列的 dense row 中 O(1) 查询。

```cpp
uint64_t sparse_dense_intersection(int i, int j) {
    if (nnz[i] > nnz[j]) std::swap(i, j);

    const auto& sp = sparse[i];
    const CountT* row = &counts[(size_t)j * KMER_DIM];

    uint64_t inter = 0;

    for (const auto& kv : sp) {
        CountT other = row[kv.id];
        inter += std::min<uint32_t>(kv.cnt, other);
    }

    return inter;
}
```

这个 kernel 避免扫描完整频次维度，适合大量短序列或稀疏序列。

### 4.4 dense SIMD kernel

若两条序列都很长、非零维度较多，直接扫描 dense row 并使用 AVX2/AVX-512 求 `min` 和累加。

AVX2 `uint16_t` 版本需要注意：`_mm256_min_epu16` 得到的是 16-bit lane，不能直接用 16-bit 累加，否则可能溢出。需要 widen 到 32-bit 后累加。

伪代码：

```cpp
uint64_t dense_simd_intersection_u16(const uint16_t* a, const uint16_t* b) {
    __m256i acc0 = _mm256_setzero_si256();
    __m256i acc1 = _mm256_setzero_si256();

    for (int d = 0; d < KMER_DIM; d += 16) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + d));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + d));
        __m256i mn = _mm256_min_epu16(va, vb);

        __m128i lo128 = _mm256_castsi256_si128(mn);
        __m128i hi128 = _mm256_extracti128_si256(mn, 1);

        __m256i lo32 = _mm256_cvtepu16_epi32(lo128);
        __m256i hi32 = _mm256_cvtepu16_epi32(hi128);

        acc0 = _mm256_add_epi32(acc0, lo32);
        acc1 = _mm256_add_epi32(acc1, hi32);
    }

    return horizontal_sum_u32(acc0) + horizontal_sum_u32(acc1);
}
```

若使用 `uint32_t` 频次，则可用 `_mm256_min_epu32`，每次处理 8 个维度。

### 4.5 kernel 自适应选择

```cpp
uint64_t intersection(int i, int j) {
    constexpr int SPARSE_THRESHOLD = 2048;

    if (std::min(nnz[i], nnz[j]) < SPARSE_THRESHOLD) {
        return sparse_dense_intersection(i, j);
    }

    return dense_simd_intersection(i, j);
}
```

`SPARSE_THRESHOLD` 需要调参，建议测试：

```text
512, 1024, 2048, 4096
```

---

## 第五阶段：Block 级并行收边与确定性 Union-Find

### 5.1 阶段目标

解决两个问题：

1. 并行计算 pair 时不能直接修改并查集，否则会有 data race；
2. 多线程发现边的顺序不确定，直接 union 会导致 parent 形态不确定，影响 MD5。

方案：并行阶段只计算和记录相似边；block 结束后，按确定的 `(i,j)` 顺序单线程执行 Union-Find。

### 5.2 为什么不用全局 local_edges + sort

朴素方案是每个线程维护 `local_edges`，最后汇总并 `std::sort`。这能保证确定性，但在高相似数据集上边数可能非常大，排序和内存占用都会成为瓶颈。

更优方案是使用 `edges_by_i`：

```cpp
struct alignas(64) EdgeRow {
    vector<int> js;
};

vector<EdgeRow> edges_by_i(N);
```

每个外层 `i` 只由一个线程处理，因此 `edges_by_i[i].js` 无需加锁。

### 5.3 Block flush 流程

```cpp
const int BLOCK_SIZE = 64;

for (int ib = 0; ib < N; ib += BLOCK_SIZE) {
    int ie = std::min(N, ib + BLOCK_SIZE);

    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = ib; i < ie; ++i) {
        auto& out = edges_by_i[i].js;
        out.clear();

        if (sum[i] == 0) continue;

        CandidateList cand = generate_candidates_adaptive(i);

        for (int j : cand) {
            if (j <= i) continue;
            if (sum[j] == 0) continue;

            if (!length_bound_may_pass(i, j, threshold)) continue;

            uint64_t inter = intersection(i, j);
            uint64_t denom = (uint64_t)sum[i] + sum[j] - inter;

            if (denom != 0 && (double)inter >= threshold * (double)denom) {
                out.push_back(j);
            }
        }

        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
    }

    for (int i = ib; i < ie; ++i) {
        for (int j : edges_by_i[i].js) {
            unite_compatible(i, j);
        }
        edges_by_i[i].js.clear();
    }
}
```

优点：

1. 不需要保存所有边；
2. 不需要全局排序；
3. union 顺序接近原始双重循环顺序；
4. 并行阶段无锁，串行 union 阶段无 data race。

### 5.4 只读 find 剪枝

在每个 block 开始时，并查集已经包含之前 block 合并过的结果。对于候选 pair，可以读取当前 root：

```cpp
if (find_no_compress(i) == find_no_compress(j)) {
    continue;
}
```

若两者已经在同一连通分量内，则无需再计算相似度，因为新的边不会改变结果。注意并行阶段不能修改 parent，只能只读判断。

---

## 第六阶段：Union-Find 兼容实现与确定性输出

### 6.1 阶段目标

保证不同线程数下输出完全一致，并通过远端脚本的 MD5 检查。

### 6.2 Union 规则复制基准程序

当前 `wj-baseline.cpp` 的合并逻辑是强制根节点为较小下标：

```cpp
if (root_i < root_j) {
    parent[root_j] = root_i;
} else {
    parent[root_i] = root_j;
}
```

优化版必须保持这一规则，不要改成 union by rank 或 union by size。否则连通分量语义可能相同，但 parent 形态不一定相同，MD5 可能失败。

推荐写法：

```cpp
int find_compress(int x) {
    if (parent[x] == x) return x;
    parent[x] = find_compress(parent[x]);
    return parent[x];
}

void unite_compatible(int a, int b) {
    int ra = find_compress(a);
    int rb = find_compress(b);

    if (ra == rb) return;

    if (ra < rb) parent[rb] = ra;
    else parent[ra] = rb;
}
```

### 6.3 全量路径压缩与输出格式

当前基准输出是一行 parent，空格分隔，最后换行：

```cpp
for (int i = 0; i < N; ++i) {
    parent[i] = find_compress(i);
}

for (int i = 0; i < N; ++i) {
    cout << parent[i] << (i == N - 1 ? "" : " ");
}
cout << '\n';
```

不要改成每行一个 parent，也不要输出序列名、调试信息或额外空白。远端脚本会直接对输出文件做 MD5。

---

## 第七阶段：编译参数、命令行与远端脚本适配

### 7.1 程序命令行兼容

程序必须支持两种运行方式：

```bash
./jaccard_cluster <fasta_file> <threshold>
./jaccard_cluster <number_of_threads> <fasta_file> <threshold>
```

远端脚本实际运行的是第二种形式：

```bash
./jaccard_cluster_test <threads> test.fasta 0.85
./jaccard_cluster_test <threads> test2.fasta 0.85
```

解析逻辑：

```cpp
if (argc == 3) {
    num_threads = 1;
    fasta_file = argv[1];
    threshold = atof(argv[2]);
} else if (argc == 4) {
    num_threads = atoi(argv[1]);
    fasta_file = argv[2];
    threshold = atof(argv[3]);
} else {
    usage();
}

omp_set_num_threads(num_threads);
```

### 7.2 Makefile 约定

当前 `Makefile` 已支持远端脚本通过 `TARGET` 覆盖输出名：

```bash
make clean TARGET=jaccard_cluster_test
make TARGET=jaccard_cluster_test
```

不要把 Makefile 改成固定只生成某个不可覆盖的二进制名。推荐继续保留：

```makefile
CXXFLAGS = -O3 -mavx2 -fopenmp -pthread
LIBS = -lz
TARGET = jaccard_cluster
SRC = wj.cpp
```

如果后续添加可选编译宏，例如 20 字母模式、dense SIMD 开关、posting 开关，也应保持远端默认命令无需额外参数即可编译通过。

### 7.3 远端执行脚本

单样例快速验证：

```powershell
.\run_test.ps1 -Runs 3 -Threads 4
```

双样例完整验证：

```powershell
.\run_all.ps1 -Runs 3 -Threads 4
```

线程数扫描：

```powershell
foreach ($t in 1,2,4,8,16,32,64) {
    .\run_all.ps1 -Runs 3 -Threads $t
}
```

脚本流程已经包含：

1. `scp` 上传 `wj.cpp` 与 `Makefile`；
2. 远端 `make clean TARGET=jaccard_cluster_test && make TARGET=jaccard_cluster_test`；
3. 远端执行 `time -p ./jaccard_cluster_test ...`；
4. 生成 `test1_out.txt` / `test2_out.txt`；
5. 与 `result1.txt` / `result2.txt` 做 `md5sum`；
6. 本地追加写入 `experiment_result.log`。

因此报告中的正确性和性能数字应优先引用 `run_all.ps1` 的输出与 `experiment_result.log`。

### 7.4 OpenMP 运行参数

如果需要测试线程绑定，可以在远端命令里加入：

```bash
export OMP_PROC_BIND=true
export OMP_PLACES=cores
```

但当前 `run_test.ps1` / `run_all.ps1` 没有暴露这两个参数。若要纳入正式实验，建议单独改脚本参数，而不是手工临时跑，避免报告数据来源混乱。

### 7.5 需要调参的参数

建议在代码顶部集中定义：

```cpp
constexpr int BLOCK_SIZE = 64;
constexpr int SPARSE_THRESHOLD = 2048;
constexpr double POSTING_ALPHA = 0.5;
constexpr int OMP_CHUNK_I = 4;
```

建议测试组合：

```text
BLOCK_SIZE: 32, 64, 128
SPARSE_THRESHOLD: 512, 1024, 2048, 4096
POSTING_ALPHA: 0.2, 0.5, 1.0
OMP_CHUNK_I: 1, 4, 8, 16
```

每次调参后至少跑：

```powershell
.\run_all.ps1 -Runs 3 -Threads 4
```

确认两个样例均 MD5 PASS 后，再做完整线程数扫描。

---

## 第八阶段：正确性验证与性能评估

### 8.1 正确性验证流程

开发过程中推荐顺序：

1. 本地只做编译级检查，确保没有明显语法问题；
2. 用 `run_test.ps1` 快速验证 `test.fasta`；
3. 用 `run_all.ps1` 验证 `test.fasta` 和 `test2.fasta`；
4. 在 `1,2,4,8,16,32,64` 线程下跑完整扫描。

正式验收命令：

```powershell
.\run_all.ps1 -Runs 3 -Threads 1
.\run_all.ps1 -Runs 3 -Threads 2
.\run_all.ps1 -Runs 3 -Threads 4
.\run_all.ps1 -Runs 3 -Threads 8
.\run_all.ps1 -Runs 3 -Threads 16
```

如果时间允许，继续跑：

```powershell
.\run_all.ps1 -Runs 3 -Threads 32
.\run_all.ps1 -Runs 3 -Threads 64
```

两个样例每次都必须 MD5 PASS。只看本地输出一致不够，最终以远端 `result1.txt` / `result2.txt` 的 MD5 对比为准。

### 8.2 若 MD5 不一致，排查顺序

1. 输出格式是否仍是一行、空格分隔、最后换行；
2. k-mer 字符集合是否和基准一致，尤其是 `X/B/Z/U/O`；
3. threshold 比较是否和基准一致；
4. union 规则是否仍是较小 root 作为根；
5. edge 收集和 union 顺序是否确定；
6. 是否出现 `uint16_t` 频次溢出；
7. SIMD 累加是否溢出或水平求和错误；
8. 并行阶段是否存在对 parent 的写操作。

### 8.3 性能数据记录

`run_all.ps1` 会输出并记录：

```text
Threads
Runs
AvgTime1
AvgTime2
MD5
```

报告中至少给出：

```text
test.fasta: 1, 2, 4, 8, 16 线程运行时间与加速比
test2.fasta: 1, 2, 4, 8, 16 线程运行时间与加速比
```

如果需要更细的性能分析，可以在远端手工使用：

```bash
perf stat ./jaccard_cluster_test 16 test.fasta 0.85 > /dev/null
```

关注：

```text
cycles
instructions
cache-misses
branch-misses
LLC-load-misses
```

但 perf / VTune 数据作为辅助分析即可，最终成绩优先看远端脚本的 MD5 与平均时间。

---

## 第九阶段：报告写法建议

报告中的优化描述建议按三个层次写。

### 9.1 算法级优化

1. 用零哈希数组替代字符串 `unordered_map`；
2. 用长度上界剪枝跳过不可能达到阈值的 pair；
3. 用 k-mer 倒排索引只生成共享至少一个 3-mer 的候选 pair；
4. 用 `sum(max)=sum_i+sum_j-sum(min)` 减少计算量。

### 9.2 数据结构与访存优化

1. dense row-major 保证每条序列的频次向量连续；
2. sparse 非零表避免扫描大量 0；
3. sparse+dense 查询和 dense SIMD 自适应选择；
4. block flush 降低边缓存内存峰值。

### 9.3 并行与确定性优化

1. k-mer 统计按序列并行，无锁；
2. pairwise 计算按外层 `i` 或 block 并行；
3. 并行阶段只记录边，不修改并查集；
4. block 结束后按确定 `(i,j)` 顺序串行 union，保证不同线程数输出一致；
5. 所有性能数据通过 `run_all.ps1` 在远端统一采集。

---

## 最终推荐优先级

### 必做优化

1. 兼容性零哈希 3-mer 统计；
2. OpenMP 并行 histogram 构建；
3. dense row-major 存储；
4. 长度上界剪枝；
5. `sum(max)=sum_i+sum_j-inter`，避免同时计算 min/max；
6. 并行只收边，串行 deterministic union；
7. 通过 `run_all.ps1` 验证两个样例 MD5。

### 强烈建议优化

1. sparse 非零特征；
2. sparse+dense intersection；
3. `edges_by_i` + block flush；
4. 倒排索引候选生成；
5. OpenMP chunk size 调参；
6. 线程数扫描并写入报告图表。

### 进阶优化

1. dense SIMD AVX2/AVX-512；
2. NUMA first-touch；
3. `OMP_PROC_BIND` 与 `OMP_PLACES`；
4. perf / VTune 指导下的软件预取。

---

## 一句话总结

本方案用兼容性零哈希替代字符串 map，用 dense+sparse 双表示降低统计和相似度计算成本，用长度上界和倒排索引减少候选 pair，用自适应 exact kernel 加速 `sum(min)`，用 block 级无锁收边和兼容 Union-Find 保证 MD5 正确性；执行与验收统一走现有 `run_test.ps1` / `run_all.ps1` 远端脚本，避免本地手工命令和正式评测流程脱节。
