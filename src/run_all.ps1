param(
    [int]$Runs = 3,
    [int]$Threads = 4
)

$Server = "root@211.87.224.231"
$Port = 8022
$Identity = "C:\Users\pc23\.ssh\id_ed25519"
$RemoteDir = "~/expr2026_baseline/src"

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "      Parallel Computing Auto-Test" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

Write-Host "`n>>> [1/4] Pushing code to server..." -ForegroundColor Yellow
scp -o StrictHostKeyChecking=no -i $Identity -P $Port wj.cpp "${Server}:${RemoteDir}/wj.cpp"
if ($LASTEXITCODE -ne 0) {
    Write-Host "Push failed for wj.cpp! Check network or ssh keys." -ForegroundColor Red
    exit 1
}

scp -o StrictHostKeyChecking=no -i $Identity -P $Port Makefile "${Server}:${RemoteDir}/Makefile"
if ($LASTEXITCODE -ne 0) {
    Write-Host "Push failed for Makefile! Check network or ssh keys." -ForegroundColor Red
    exit 1
}

Write-Host ">>> [2/4] Compiling on server..." -ForegroundColor Yellow
$compileCmd = "cd $RemoteDir && make clean TARGET=jaccard_cluster_test && make TARGET=jaccard_cluster_test"
ssh -o StrictHostKeyChecking=no -i $Identity -p $Port $Server $compileCmd

if ($LASTEXITCODE -ne 0) {
    Write-Host "Compilation failed! Check your code syntax." -ForegroundColor Red
    exit 1
}

$totalRealTime1 = 0.0
$totalRealTime2 = 0.0

Write-Host ">>> [3/4] Running tests (Total $Runs runs)..." -ForegroundColor Yellow
for ($i = 1; $i -le $Runs; $i++) {
    Write-Host "`n  --- Run $i ---"
    
    # --- TEST 1 ---
    Write-Host "  Testing Sample 1 (test.fasta, threshold 0.85) ..."
    $runCmd1 = "cd $RemoteDir && { time -p ./jaccard_cluster_test $Threads test.fasta 0.85 > test1_out.txt; } 2>&1 && md5sum test1_out.txt && md5sum result1.txt"
    $output1 = ssh -o StrictHostKeyChecking=no -i $Identity -p $Port $Server $runCmd1
    
    $testMd5_1 = $null
    $baseMd5_1 = $null
    foreach ($line in $output1) {
        if ($line -match "test1_out.txt") { $testMd5_1 = ($line -split '\s+')[0] }
        if ($line -match "result1.txt") { $baseMd5_1 = ($line -split '\s+')[0] }
    }
    
    if ($testMd5_1 -ne $baseMd5_1) {
        Write-Host "    [ERROR] Sample 1 MD5 mismatch! Output is incorrect." -ForegroundColor Red
        exit 1
    }
    $realTimeLine1 = $output1 | Select-String -Pattern "^real\s+"
    if ($realTimeLine1) {
        $realTime1 = [double]($realTimeLine1.ToString() -split '\s+')[1]
        Write-Host "    -> Sample 1 Time: $realTime1 seconds, MD5: PASS" -ForegroundColor Green
        $totalRealTime1 += $realTime1
    }

    # --- TEST 2 ---
    Write-Host "  Testing Sample 2 (test2.fasta, threshold 0.85) ..."
    $runCmd2 = "cd $RemoteDir && { time -p ./jaccard_cluster_test $Threads test2.fasta 0.85 > test2_out.txt; } 2>&1 && md5sum test2_out.txt && md5sum result2.txt"
    $output2 = ssh -o StrictHostKeyChecking=no -i $Identity -p $Port $Server $runCmd2
    
    $testMd5_2 = $null
    $baseMd5_2 = $null
    foreach ($line in $output2) {
        if ($line -match "test2_out.txt") { $testMd5_2 = ($line -split '\s+')[0] }
        if ($line -match "result2.txt") { $baseMd5_2 = ($line -split '\s+')[0] }
    }
    
    if ($testMd5_2 -ne $baseMd5_2) {
        Write-Host "    [ERROR] Sample 2 MD5 mismatch! Output is incorrect." -ForegroundColor Red
        Write-Host "    Expected : $baseMd5_2"
        Write-Host "    Actual   : $testMd5_2"
        exit 1
    }
    $realTimeLine2 = $output2 | Select-String -Pattern "^real\s+"
    if ($realTimeLine2) {
        $realTime2 = [double]($realTimeLine2.ToString() -split '\s+')[1]
        Write-Host "    -> Sample 2 Time: $realTime2 seconds, MD5: PASS" -ForegroundColor Green
        $totalRealTime2 += $realTime2
    }
}

$avgTime1 = [math]::Round($totalRealTime1 / $Runs, 3)
$avgTime2 = [math]::Round($totalRealTime2 / $Runs, 3)

Write-Host "`n>>> [4/4] Test Report" -ForegroundColor Cyan
Write-Host "  ======================================"
Write-Host "  All $Runs tests passed MD5 check for BOTH samples." -ForegroundColor Green
Write-Host "  OMP_NUM_THREADS : $Threads"
Write-Host "  Sample 1 Avg Time: $avgTime1 seconds" -ForegroundColor Green
Write-Host "  Sample 2 Avg Time: $avgTime2 seconds" -ForegroundColor Green
Write-Host "  ======================================"

$logLine = "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] Threads: $Threads | Runs: $Runs | AvgTime1: ${avgTime1}s | AvgTime2: ${avgTime2}s | MD5: PASS"
$logLine | Out-File -FilePath "experiment_result.log" -Append
Write-Host "Saved to experiment_result.log"
