$exe = "build\source\Release\ParticleSimulation.exe"
$scriptDir = "script\"
$outputDir = "output"

if (-not (Test-Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir | Out-Null
}

$simulations = @(
    @{ arg = "anisotropy"; script = "anisotropy.py" },
    @{ arg = "finite-diff-err"; script = "finite_diff_rel_error.py" },
    @{ arg = "assignment-schemes"; script = "assignment_schemes.py" },
    @{ arg = "poor-man-vs-laplacian"; script = "poor_man_vs_laplace.py" },
    @{ arg = "pm-accuracy-soft"; script = "pm_accuracy.py"; extraArgs = @("soft") },
    @{ arg = "pm-accuracy-no-soft"; script = "pm_accuracy.py"; extraArgs = @("no-soft") },
    @{ arg = "pm-optimal"; script = "pm_optimal.py" },
    @{ arg = "pm-optimal-var-diameter"; script = "pm_ref_approx_qual.py" },
    @{ arg = "pm-optimal-assignments"; script = "optimal_assignments.py" },
    @{ arg = "p3m-accuracy-assignments"; script = "p3m_accuracy.py" },
    @{ arg = "p3m-accuracy-shapes"; script = "p3m_accuracy_shapes.py" },
    @{ arg = "bh-accuracy"; script = "bh_accuracy.py" }
)

foreach ($sim in $simulations) {
    Write-Host "Running simulation: $($sim.arg)"

    & $exe $sim.arg $outputDir *> $null

    $scriptPath = Join-Path $scriptDir $sim.script

    if ($sim.ContainsKey("extraArgs")) {
        & python $scriptPath @($sim.extraArgs) *> $null
    } else {
        & python $scriptPath *> $null
    }
}
