param (
    [Parameter(Mandatory = $true)]
    [string]$SimulationName
)

$exe = "build\source\Release\ParticleSimulation.exe"
$scriptDir = "script\"
$animateScript = Join-Path $scriptDir "animate2d.py"
$conservedQuantitiesScript = Join-Path $scriptDir "conserved_quantities.py"
$outputDir = "output"

if (-not (Test-Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir | Out-Null
}

$simulations = @(
    @{ arg = "galaxy-sim-pm"; objectsCnt = "1" },
    @{ arg = "galaxy-sim-p3m"; objectsCnt = "1" },
    @{ arg = "galaxy-sim-bh"; objectsCnt = "1" },
    @{ arg = "cluster-sim-bh"; objectsCnt = "1" },
    @{ arg = "cluster-sim-pm"; objectsCnt = "1" },
    @{ arg = "cluster-sim-p3m"; objectsCnt = "1" },
    @{ arg = "galaxy-collision-sim-bh"; objectsCnt = "2" },
    @{ arg = "galaxy-collision-sim-pm"; objectsCnt = "2" },
    @{ arg = "galaxy-collision-sim-p3m"; objectsCnt = "2" }
)

$selectedSim = $simulations | Where-Object { $_.arg -eq $SimulationName }

if (-not $selectedSim) {
    Write-Host "Invalid simulation name: $SimulationName"
    Write-Host "Valid options are:"
    $simulations | ForEach-Object { Write-Host "  - $($_.arg)" }
    exit 1
}

# Run simulation
& $exe $selectedSim.arg $outputDir
if ($LASTEXITCODE -ne 0) {
    Write-Error "Executable failed with exit code $LASTEXITCODE. Terminating script."
    exit $LASTEXITCODE
}

# Run Python scripts
& python $animateScript $selectedSim.objectsCnt
& python $conservedQuantitiesScript
