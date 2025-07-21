@echo off
cd /d "%~dp0"
Rshiny\R-Portable\App\R-Portable\bin\Rscript.exe -e "shiny::runApp('Rshiny/app.R', launch.browser=TRUE)"