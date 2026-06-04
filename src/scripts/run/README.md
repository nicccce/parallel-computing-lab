# 运行脚本说明

本目录保留两个 PowerShell 运行入口：`run_test.ps1` 和 `run_all.ps1`。脚本会自动定位项目根目录，所以可以从项目根目录运行，也可以进入 `scripts/run` 后运行。

## 前置条件

- 本机可以通过 `ssh` / `scp` 访问远端服务器。
- SSH 私钥路径为 `C:\Users\pc23\.ssh\id_ed25519`。
- 远端实验目录为 `~/expr2026_baseline/src`，其中应已有 `test.fasta`、`test2.fasta`、`result1.txt`。
- 远端还应已有 `test2.fasta` 多阈值参考答案，默认路径为 `~/expr2026_baseline/src/baseline_answers/test2_threshold_*.txt`。
- 远端参考答案必须使用 Linux LF 换行；如果从 Windows 上传了 CRLF 文件，即使内容一致也会导致 raw MD5 mismatch。
- 本地根目录保留 `wj.cpp` 和 `Makefile`。

## run_test.ps1

用途：快速验证 `test.fasta` 单个数据集。

推荐从项目根目录运行：

```powershell
.\scripts\run\run_test.ps1 -Runs 3 -Threads 4
```

脚本会做这些事：

1. 上传本地 `wj.cpp` 和 `Makefile` 到远端实验目录。
2. 在远端执行 `make clean TARGET=jaccard_cluster_test` 和 `make TARGET=jaccard_cluster_test`。
3. 远端重复运行 `./jaccard_cluster_test <Threads> test.fasta 0.85`。
4. 将程序输出和远端 `result1.txt` 做 `md5sum` 对比。
5. 将平均运行时间和 MD5 状态追加到 `results/experiment_result.log`。

## run_all.ps1

用途：完整验证 `test.fasta` 和 `test2.fasta`。

推荐从项目根目录运行：

```powershell
.\scripts\run\run_all.ps1 -Runs 3 -Threads 4
```

`test2.fasta` 的默认阈值为 `0.80, 0.85, 0.90, 0.95`。也可以指定：

```powershell
.\scripts\run\run_all.ps1 -Runs 3 -Threads 4 -Test2Thresholds 0.80,0.85,0.90,0.95
```

脚本会做这些事：

1. 只上传本地 `wj.cpp` 和 `Makefile` 到远端实验目录，不上传答案文件。
2. 在远端重新编译 `jaccard_cluster_test`。
3. 对 `test.fasta` 固定阈值 `0.85`，重复运行 `-Runs` 次，并与远端 `result1.txt` 做 MD5 校验，最后计算平均时间。
4. 对 `test2.fasta` 的每个阈值各运行一次，并与远端已有的 `baseline_answers/test2_threshold_*.txt` 做 MD5 校验。
5. 对 `test2.fasta` 的多个阈值运行时间求平均。
6. 将结果追加到 `results/experiment_result.log`。

线程数扫描示例：

```powershell
foreach ($t in 1,2,4,8,16,32,64) {
    .\scripts\run\run_all.ps1 -Runs 3 -Threads $t
}
```

## 参数

- `-Runs`：`test.fasta` 重复运行次数，默认 `3`。`test2.fasta` 不按 Runs 重复，而是每个阈值运行一次。
- `-Threads`：传给程序的线程数，默认 `4`。
- `-Test2Thresholds`：`test2.fasta` 的阈值列表，默认 `0.80, 0.85, 0.90, 0.95`。
- `-RemoteTest2AnswerDir`：远端 `test2.fasta` 多阈值答案目录，默认 `baseline_answers`，相对于 `~/expr2026_baseline/src`。

## 输出位置

- 本地运行日志：`results/experiment_result.log`
- 远端参考答案：`baseline_answers/test2_threshold_*.txt`
- 远端临时输出：`test1_out.txt`、`test2_out_<threshold>.txt`
- 远端可执行文件：`jaccard_cluster_test`

## 常用命令

如果 Windows PowerShell 提示“在此系统上禁止运行脚本”，不要改全局执行策略，直接用一次性的 `-ExecutionPolicy Bypass -File` 启动即可：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\run\run_test.ps1 -Runs 1 -Threads 4
powershell -ExecutionPolicy Bypass -File .\scripts\run\run_all.ps1 -Runs 1 -Threads 4
```

建议每次改完 `wj.cpp` 后先做一次 1 轮冒烟测试：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\run\run_test.ps1 -Runs 1 -Threads 4
```

确认 `test.fasta` 的 `MD5: PASS` 后，再跑完整验证：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\run\run_all.ps1 -Runs 1 -Threads 4
```

若脚本报 MD5 mismatch，优先检查 `wj.cpp` 输出格式、并查集合并规则、阈值判断和 k-mer 字符映射是否与基准程序一致。
