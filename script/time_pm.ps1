$exePath = "C:/Projects/ParticleSimulation/build/source/Release/ParticleSimulation.exe"

$outputFile = "C:/Projects/ParticleSimulation/pm-timing.txt"

"Parameter N, Execution Time (ms)" | Out-File -FilePath $outputFile -Encoding UTF8

for ($N = 50000; $N -le 1000000; $N += 50000) {
    Write-Host "Running with N = $N"

    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()

    & $exePath $N *> $null

    $stopwatch.Stop()
    $elapsedMs = $stopwatch.Elapsed.TotalMilliseconds

    "$N, $elapsedMs" | Out-File -FilePath $outputFile -Encoding UTF8 -Append
}
