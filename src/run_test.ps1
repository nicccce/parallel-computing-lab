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

$totalRealTime = 0.0

Write-Host ">>> [3/4] Running tests (Total $Runs runs)..." -ForegroundColor Yellow
for ($i = 1; $i -le $Runs; $i++) {
    Write-Host "  Running test $i ..."
    
    $runCmd = "cd $RemoteDir && { time -p ./jaccard_cluster_test $Threads test.fasta 0.85 > test_output.txt; } 2>&1 && md5sum test_output.txt && md5sum result1.txt"
    $output = ssh -o StrictHostKeyChecking=no -i $Identity -p $Port $Server $runCmd
    
    $testMd5 = $null
    $baseMd5 = $null
    foreach ($line in $output) {
        if ($line -match "test_output.txt") { $testMd5 = ($line -split '\s+')[0] }
        if ($line -match "result1.txt") { $baseMd5 = ($line -split '\s+')[0] }
    }
    
    if ($testMd5 -ne $baseMd5) {
        Write-Host "    [ERROR] MD5 mismatch! Output is incorrect." -ForegroundColor Red
        Write-Host "    Expected : $baseMd5"
        Write-Host "    Actual   : $testMd5"
        exit 1
    }
    
    $realTimeLine = $output | Select-String -Pattern "^real\s+"
    if ($realTimeLine) {
        $realTimeStr = ($realTimeLine.ToString() -split '\s+')[1]
        $realTime = [double]$realTimeStr
        Write-Host "    -> Time: $realTime seconds, MD5: PASS" -ForegroundColor Green
        $totalRealTime += $realTime
    } else {
        Write-Host "    [WARNING] Could not parse real time!" -ForegroundColor Yellow
        $output
    }
}

$avgTime = $totalRealTime / $Runs
$avgTimeStr = [math]::Round($avgTime, 3)

Write-Host "`n>>> [4/4] Test Report" -ForegroundColor Cyan
Write-Host "  ======================================"
Write-Host "  All $Runs tests passed MD5 check." -ForegroundColor Green
Write-Host "  OMP_NUM_THREADS : $Threads"
Write-Host "  Avg Real Time   : $avgTimeStr seconds" -ForegroundColor Green
Write-Host "  ======================================"

$logLine = "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] Threads: $Threads | Runs: $Runs | AvgTime: ${avgTimeStr}s | MD5: PASS"
$logLine | Out-File -FilePath "experiment_result.log" -Append
Write-Host "Saved to experiment_result.log"
