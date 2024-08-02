#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
#setenv SCRAM_ARCH slc6_amd64_gcc700
export SCRAM_ARCH=slc7_amd64_gcc10
scramv1 project CMSSW CMSSW_12_4_13
cd CMSSW_12_4_13/src
eval `scramv1 runtime -sh`
cd ${_CONDOR_SCRATCH_DIR}/CMSSW_12_4_13/src/
cp ${_CONDOR_SCRATCH_DIR}/*.py .
cp ${_CONDOR_SCRATCH_DIR}/*.json .
export X509_USER_PROXY=/afs/cern.ch/user/b/bgreenbe/x509up_u104084
#scram b -j 4
#eval `scramv1 runtime -csh`
mv ../../flist_*.txt .
#ls
#FOR DATA
python3 plot_2mu2e.py $2 $3 $4
#python3 plot_2mu2e_backup.py $2 $3 $4
echo "Done with args $2 $3 $4."
