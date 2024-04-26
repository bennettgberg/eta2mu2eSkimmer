#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc10
scramv1 project CMSSW CMSSW_12_4_13
cd CMSSW_12_4_13/src
eval `scramv1 runtime -sh`
cd ${_CONDOR_SCRATCH_DIR}/CMSSW_12_4_13/src/
cp ${_CONDOR_SCRATCH_DIR}/*.py .
cp ${_CONDOR_SCRATCH_DIR}/*.C .
export X509_USER_PROXY=/afs/cern.ch/user/b/bgreenbe/x509up_u104084
#scram b -j 4
#eval `scramv1 runtime -csh`
mkdir trigfilelists
mv ../../*.txt ./trigfilelists/
#ls
#echo "trigfilelists:"
#ls trigfilelists
#FOR DATA
echo root -q -l -b -e ".L measure_trigCorrex.C+" -e "measure_trigCorrex(\"$2\", \"$3\", \"$4\")"
root -q -l -b -e ".L measure_trigCorrex.C+" -e "measure_trigCorrex(\"$2\", \"$3\", \"$4\")"
echo "Done with args $2 $3 $4."
