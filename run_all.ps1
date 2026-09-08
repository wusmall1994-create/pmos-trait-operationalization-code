$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location -LiteralPath $projectRoot
try {
    $rscript = Get-Command Rscript -ErrorAction Stop
    & $rscript.Source --vanilla (Join-Path $projectRoot 'run_all.R')
    if ($LASTEXITCODE -ne 0) {
        throw "R pipeline exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

