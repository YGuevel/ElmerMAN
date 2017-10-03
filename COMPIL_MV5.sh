#!/bin/sh
# YG@2016
# yann.guevel@univ-ubs.fr
# Part of ANM continuation/detection/branch-switching in ELMER FEM
# Jean-Marc Cadou / Gregory girault / Yann Guevel
# Institut de Recherche Dupuy de Lome - Universite Bretagne Sud



# $1 == nom du fichier solver file.f90
filetocompile=$1


# ENVIRONNEMENT ELMER avec MUMPS V5 + OpenBLAS
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/lab/limatb/app/MUMPS/MUMPS_5.0.0_FILES_fPIC_OpenBLAS/lib
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/lab/limatb/app/MUMPS/OpenBLAS/Install/lib
#setenv ELMER /lab/limatb/app/ELMER/VERSIONs/ELMERDEB7MUMPSV5OpenBLAS
#setenv PATH $ELMER/bin:${PATH}
#setenv OPENBLAS_NUM_THREADS 16


module1=MOD-UtilsFEM
module2=MOD-HomeMadeDirectSolvers-v01
module3=MOD-DiscreteOperators-v03
module4=MOD-ELMERModResOut-v03
module5=MOD-ANMToolBox-V05
module6=MOD-FlowsolveCorrectionClassic  
module7=MOD-ANM-Singularities


MODULESf90=$module1'.f90 '$module2'.f90 '$module3'.f90 '$module4'.f90 '$module5'.f90 '$module6'.f90 '$module7'.f90'
MODULESSO=$module1'.o '$module2'.o '$module3'.o '$module4'.o '$module5'.o '$module6'.o '$module7'.o'


# LIBMPISEQ="/lab/home/limatb/guevel/_HOME_PRO_/ELMER/MUMPS/MUMPS_5.0.0/libseq"
# LIBMUMPS="/lab/home/limatb/guevel/_HOME_PRO_/ELMER/MUMPS/MUMPS_5.0.0/libseq"
# OPTIONS="-L$LIBMPISEQ -lmpiseq -L$LIBMUMPS -ldmumps -lmumps_common -lpord"

#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/lab/limatb/app/MUMPS/MUMPS_5.0.0_FILES_fPIC_OpenBLAS/lib
#setenv ELMER /lab/limatb/app/ELMER/VERSIONs/ELMERDEB7MUMPSV5OpenBLAS
#setenv PATH $ELMER/bin:${PATH}

# setenv OPENBLAS_NUM_THREADS 16
# By default, if there are 2 hyper-threading logical cores, OpenBLAS will use only one core.
# You can compile OpenBLAS with NO_AFFINITY=1 to disable this feature

# elmerf90 -c $MODULESf90 $OPTIONS -fPIC
ls $MODULESf90
elmerf90 -c $MODULESf90
ls $MODULESSO
# elmerf90 $MODULESSO $filetocompile $OPTIONS -fPIC -o ${filetocompile%.f90}.so
elmerf90 $MODULESSO $filetocompile -o ${filetocompile%.f90}Mv5.so
