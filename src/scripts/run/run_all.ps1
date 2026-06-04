param(
    [int]$Runs = 3,
    [int]$Threads = 4,
    [string[]]$Test2Thresholds = @("0.80", "0.85", "0.90", "0.95"),
    [string]$RemoteTest2AnswerDir = "baseline_answers"
)

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..\..")).Path
$ResultsDir = Join-Path $RepoRoot "results"
$SourceFile = Join-Path $RepoRoot "wj.cpp"
$MakefilePath = Join-Path $RepoRoot "Makefile"
$LogPath = Join-Path $ResultsDir "experiment_result.log"

New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null

$Test2Thresholds = @(
    $Test2Thresholds |
        ForEach-Object { $_ -split "," } |
        ForEach-Object { $_.Trim() } |
        Where-Object { $_ -ne "" }
)

if ($Test2Thresholds.Count -eq 0) {
    Write-Host "No test2 thresholds provided." -ForegroundColor Red
    exit 1
}

if (-not (Test-Path -LiteralPath $SourceFile)) {
    Write-Host "Cannot find wj.cpp at $SourceFile" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path -LiteralPath $MakefilePath)) {
    Write-Host "Cannot find Makefile at $MakefilePath" -ForegroundColor Red
    exit 1
}

function Get-ThresholdTag([string]$Threshold) {
    return ($Threshold -replace "\.", "_")
}

function Get-Test2ReferenceName([string]$Threshold) {
    return "test2_threshold_$(Get-ThresholdTag $Threshold).txt"
}

$Server = "root@211.87.224.231"
$Port = 8022
$Identity = "C:\Users\pc23\.ssh\id_ed25519"
$RemoteDir = "~/expr2026_baseline/src"
$SshOptions = @(
    "-o", "StrictHostKeyChecking=no",
    "-o", "BatchMode=yes",
    "-o", "ConnectTimeout=20",
    "-o", "ServerAliveInterval=10",
    "-o", "ServerAliveCountMax=3"
)
$ScpOptions = @(
    "-o", "StrictHostKeyChecking=no",
    "-o", "BatchMode=yes",
    "-o", "ConnectTimeout=20",
    "-o", "ServerAliveInterval=10",
    "-o", "ServerAliveCountMax=3"
)

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "      Parallel Computing Auto-Test" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Threads          : $Threads"
Write-Host "test.fasta Runs  : $Runs"
Write-Host "test2 Thresholds : $($Test2Thresholds -join ', ')"
Write-Host "test2 Answers    : $RemoteTest2AnswerDir/test2_threshold_<threshold>.txt on server"

Write-Host "`n>>> [1/5] Pushing code to server..." -ForegroundColor Yellow
scp @ScpOptions -i $Identity -P $Port $SourceFile "${Server}:${RemoteDir}/wj.cpp"
if ($LASTEXITCODE -ne 0) {
    Write-Host "Push failed for wj.cpp! Check network or ssh keys." -ForegroundColor Red
    exit 1
}

scp @ScpOptions -i $Identity -P $Port $MakefilePath "${Server}:${RemoteDir}/Makefile"
if ($LASTEXITCODE -ne 0) {
    Write-Host "Push failed for Makefile! Check network or ssh keys." -ForegroundColor Red
    exit 1
}

Write-Host ">>> [2/5] Compiling on server..." -ForegroundColor Yellow
$compileCmd = "cd $RemoteDir && make clean TARGET=jaccard_cluster_test && make TARGET=jaccard_cluster_test"
ssh @SshOptions -i $Identity -p $Port $Server $compileCmd

if ($LASTEXITCODE -ne 0) {
    Write-Host "Compilation failed! Check your code syntax." -ForegroundColor Red
    exit 1
}

$totalRealTime1 = 0.0
$totalRealTime2 = 0.0
$test2Rows = @()

Write-Host ">>> [3/5] Running test.fasta (threshold 0.85, total $Runs runs)..." -ForegroundColor Yellow
for ($i = 1; $i -le $Runs; $i++) {
    Write-Host "  Run $i ..."

    $runCmd1 = "cd $RemoteDir && { time -p ./jaccard_cluster_test $Threads test.fasta 0.85 > test1_out.txt; } 2>&1 && md5sum test1_out.txt && md5sum result1.txt"
    $output1 = ssh @SshOptions -i $Identity -p $Port $Server $runCmd1

    $testMd5_1 = $null
    $baseMd5_1 = $null
    foreach ($line in $output1) {
        if ($line -match "test1_out.txt") { $testMd5_1 = ($line -split '\s+')[0] }
        if ($line -match "result1.txt") { $baseMd5_1 = ($line -split '\s+')[0] }
    }

    if ($testMd5_1 -ne $baseMd5_1) {
        Write-Host "    [ERROR] test.fasta MD5 mismatch! Output is incorrect." -ForegroundColor Red
        Write-Host "    Expected : $baseMd5_1"
        Write-Host "    Actual   : $testMd5_1"
        exit 1
    }

    $realTimeLine1 = $output1 | Select-String -Pattern "^real\s+"
    if (-not $realTimeLine1) {
        Write-Host "    [ERROR] Could not parse test.fasta real time." -ForegroundColor Red
        $output1
        exit 1
    }

    $realTime1 = [double]($realTimeLine1.ToString() -split '\s+')[1]
    Write-Host "    -> test.fasta Time: $realTime1 seconds, MD5: PASS" -ForegroundColor Green
    $totalRealTime1 += $realTime1
}

Write-Host ">>> [4/5] Running test2.fasta thresholds once..." -ForegroundColor Yellow
foreach ($threshold in $Test2Thresholds) {
    $tag = Get-ThresholdTag $threshold
    $outName = "test2_out_${tag}.txt"
    $refName = Get-Test2ReferenceName $threshold

    Write-Host "  Threshold $threshold ..."
    $remoteRefPath = "$RemoteTest2AnswerDir/$refName"
    $runCmd2 = "cd $RemoteDir && { time -p ./jaccard_cluster_test $Threads test2.fasta $threshold > $outName; } 2>&1 && md5sum $outName && md5sum $remoteRefPath"
    $output2 = ssh @SshOptions -i $Identity -p $Port $Server $runCmd2

    $testMd5_2 = $null
    $baseMd5_2 = $null
    foreach ($line in $output2) {
        if ($line -match [regex]::Escape($outName)) { $testMd5_2 = ($line -split '\s+')[0] }
        if ($line -match [regex]::Escape($remoteRefPath)) { $baseMd5_2 = ($line -split '\s+')[0] }
    }

    if ($testMd5_2 -ne $baseMd5_2) {
        Write-Host "    [ERROR] test2.fasta threshold $threshold MD5 mismatch!" -ForegroundColor Red
        Write-Host "    Reference: $remoteRefPath"
        Write-Host "    Expected : $baseMd5_2"
        Write-Host "    Actual   : $testMd5_2"
        exit 1
    }

    $realTimeLine2 = $output2 | Select-String -Pattern "^real\s+"
    if (-not $realTimeLine2) {
        Write-Host "    [ERROR] Could not parse test2.fasta threshold $threshold real time." -ForegroundColor Red
        $output2
        exit 1
    }

    $realTime2 = [double]($realTimeLine2.ToString() -split '\s+')[1]
    Write-Host "    -> test2.fasta threshold $threshold Time: $realTime2 seconds, MD5: PASS" -ForegroundColor Green
    $totalRealTime2 += $realTime2
    $test2Rows += [pscustomobject]@{
        Threshold = $threshold
        Seconds = $realTime2
    }
}

$avgTime1 = [math]::Round($totalRealTime1 / $Runs, 3)
$avgTime2 = [math]::Round($totalRealTime2 / $Test2Thresholds.Count, 3)

Write-Host "`n>>> [5/5] Test Report" -ForegroundColor Cyan
Write-Host "  ======================================"
Write-Host "  test.fasta passed $Runs MD5 checks at threshold 0.85." -ForegroundColor Green
Write-Host "  test2.fasta passed all threshold MD5 checks." -ForegroundColor Green
Write-Host "  OMP_NUM_THREADS : $Threads"
Write-Host "  test.fasta Avg Time: $avgTime1 seconds" -ForegroundColor Green
Write-Host "  test2.fasta Threshold Times:"
$test2Rows | Format-Table -AutoSize
Write-Host "  test2.fasta Avg Time Across Thresholds: $avgTime2 seconds" -ForegroundColor Green
Write-Host "  ======================================"

$thresholdText = $Test2Thresholds -join "/"
$logLine = "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] Threads: $Threads | test.fasta Runs: $Runs | test.fasta AvgTime: ${avgTime1}s | test2 Thresholds: $thresholdText | test2 AvgTimeAcrossThresholds: ${avgTime2}s | MD5: PASS"
$logLine | Out-File -FilePath $LogPath -Append -Encoding Unicode
Write-Host "Saved to $LogPath"
