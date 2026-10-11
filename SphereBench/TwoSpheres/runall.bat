:: run all 3 problems on all 4 boundary types
:: problem 3 has 3 separate cases, two of which require quarters=2
:: assumes Gmsh and GetDP are in the path
:: script with coarse mesh (s=1.5), can be changed below

@echo off

set "OPT=-setnumber s 1.5 -setnumber quarters 1 -setnumber mags 1"
gmsh main.geo -3 %OPT%
call :getdp 1
call :getdp 2
call :getdp 3
set "OPT=-setnumber s 1.5 -setnumber quarters 2 -setnumber mags 2"
gmsh main.geo -3 %OPT%
call :getdp 3
set "OPT=-setnumber s 1.5 -setnumber quarters 2 -setnumber mags 3"
call :getdp 3
goto :eof

:getdp
set "PROB=-setnumber prob %~1"
for /l %%b in (1,1,4) do (
    getdp main.pro -solve ResMain -pos PostMain %OPT% %PROB% -setnumber bound %%b
)
exit /B

:eof