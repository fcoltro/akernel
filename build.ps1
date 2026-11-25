# Configure, build and test akernel with the MSVC Build Tools.
#
#   .\build.ps1                 Release build, run all tests
#   .\build.ps1 -Config Debug   Debug build, run all tests
#   .\build.ps1 -Bench          also run the benchmark program
#   .\build.ps1 -Filter ssi     run only tests whose name matches
#   .\build.ps1 -Workbench      build and start the desktop workbench
#                               (apps/workbench; fetches GLFW and Dear ImGui
#                               once, needs git), no tests
param(
    [ValidateSet('Release', 'Debug')] [string]$Config = 'Release',
    [switch]$Bench,
    [string]$Filter = '',
    [switch]$Workbench
)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot

# Locate vcvars64.bat through vswhere and import its environment once.
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path $vswhere)) { throw 'vswhere.exe not found; install the Visual Studio Build Tools.' }
$vsPath = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
$vcvars = Join-Path $vsPath 'VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path $vcvars)) { throw "vcvars64.bat not found under $vsPath" }

cmd /c "`"$vcvars`" >nul 2>nul && set" | ForEach-Object {
    if ($_ -match '^([^=]+)=(.*)$') { Set-Item -Path "env:$($Matches[1])" -Value $Matches[2] }
}

$buildDir = Join-Path $root "build\$Config"
$workbenchFlag = if ($Workbench) { 'ON' } else { 'OFF' }
cmake -S $root -B $buildDir -G Ninja "-DCMAKE_BUILD_TYPE=$Config" "-DAK_BUILD_WORKBENCH=$workbenchFlag" | Out-Host
if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }

if ($Workbench) {
    cmake --build $buildDir --target ak_workbench | Out-Host
    if ($LASTEXITCODE -ne 0) { throw 'Workbench build failed.' }
    $exe = Join-Path $buildDir 'apps\workbench\ak_workbench.exe'
    Write-Host "Starting $exe"
    Start-Process -FilePath $exe -WorkingDirectory $root
    exit 0
}
cmake --build $buildDir | Out-Host
if ($LASTEXITCODE -ne 0) { throw 'Build failed.' }

$ctestArgs = @('--test-dir', $buildDir, '--output-on-failure')
if ($Filter) { $ctestArgs += @('-R', $Filter) }
ctest @ctestArgs | Out-Host
$testExit = $LASTEXITCODE

if ($Bench) {
    $benchExe = Join-Path $buildDir 'ak_bench.exe'
    if (Test-Path $benchExe) { & $benchExe | Out-Host } else { Write-Host 'ak_bench.exe was not built.' }
}

exit $testExit
