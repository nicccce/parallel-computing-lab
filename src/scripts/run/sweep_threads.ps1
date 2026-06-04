param(
    [string]$ThreadList = "1,2,8,16,32,64",
    [int]$Runs = 1
)

$threads = $ThreadList -split "," | ForEach-Object { [int]$_.Trim() }

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

foreach ($t in $threads) {
    Write-Host ""
    Write-Host "=========================================" -ForegroundColor Cyan
    Write-Host "  Thread Sweep: $t threads" -ForegroundColor Cyan
    Write-Host "=========================================" -ForegroundColor Cyan

    & "$ScriptDir\run_all.ps1" -Runs $Runs -Threads $t

    if ($LASTEXITCODE -ne 0) {
        Write-Host "FAILED at $t threads, stopping sweep." -ForegroundColor Red
        exit 1
    }
}

Write-Host ""
Write-Host "All thread counts completed." -ForegroundColor Green
