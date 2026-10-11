:: run all 7 problems on all 2 boundary types
:: assumes Gmsh and GetDP are in the path
:: script with coarse mesh (s=1.5), can be changed below

@echo off

set "OPT=-setnumber s 1.5"
gmsh main.geo -3 %OPT%
for /l %%p in (1,1,7) do (
    call :getdp %%p
)
goto :eof

:getdp
set "PROB=-setnumber prob %~1"
getdp main.pro -solve ResMain -pos PostMain %OPT% %PROB% -setnumber iabc 1
getdp main.pro -solve ResMain -pos PostMain %OPT% %PROB% -setnumber iabc 0
exit /B

:eof