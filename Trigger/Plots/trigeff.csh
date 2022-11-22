#!/bin/csh

#pyscript=$1
#inputfile=$2
#outputname=$3
#Jet/Displaced Trigger = $4
#Data/MC=$5
#Control/Signal=$6
#Blind/Not=$7

source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv ROOTSYS /cvmfs/cms.cern.ch/slc6_amd64_gcc472/lcg/root/5.32.00-cms
setenv PATH $ROOTSYS/bin:$PATH
set current=`pwd`
echo $current
scram p CMSSW CMSSW_7_6_2
cd CMSSW_7_6_2
cmsenv
cd $current

setenv LD_LIBRARY_PATH $ROOTSYS/lib/:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
python $1 -f $2 -o $3 -w $4 $5 
