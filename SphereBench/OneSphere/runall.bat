# run all 8 problems on all 4 boundary types
# assumes Gmsh and GetDP are in the path
# script with coarse mesh (s=1.5), can be changed below

set "OPT=-setnumber s 1.5"
gmsh main.geo -3 %OPT%
for /l %%p in (1,1,8) do (
    call :getdp %%p
)
goto :eof

:getdp
set "PROB=-setnumber prob %~1"
for /l %%b in (1,1,4) do (
    getdp main.pro -solve ResMain -pos PostMain %OPT% %PROB% -setnumber bound %%b
)
exit /B

:eof