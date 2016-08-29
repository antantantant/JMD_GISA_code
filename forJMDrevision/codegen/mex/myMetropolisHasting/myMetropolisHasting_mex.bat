@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2013a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2013a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=myMetropolisHasting_mex
set MEX_NAME=myMetropolisHasting_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for myMetropolisHasting > myMetropolisHasting_mex.mki
echo COMPILER=%COMPILER%>> myMetropolisHasting_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> myMetropolisHasting_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> myMetropolisHasting_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> myMetropolisHasting_mex.mki
echo LINKER=%LINKER%>> myMetropolisHasting_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> myMetropolisHasting_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> myMetropolisHasting_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> myMetropolisHasting_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> myMetropolisHasting_mex.mki
echo BORLAND=%BORLAND%>> myMetropolisHasting_mex.mki
echo OMPFLAGS= >> myMetropolisHasting_mex.mki
echo OMPLINKFLAGS= >> myMetropolisHasting_mex.mki
echo EMC_COMPILER=msvcsdk>> myMetropolisHasting_mex.mki
echo EMC_CONFIG=optim>> myMetropolisHasting_mex.mki
"C:\Program Files\MATLAB\R2013a\bin\win64\gmake" -B -f myMetropolisHasting_mex.mk
