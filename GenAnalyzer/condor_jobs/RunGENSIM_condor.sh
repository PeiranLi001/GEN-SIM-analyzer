#!/bin/bash

echo "Starting job on " `date` #Only to display the starting of production date
echo "Running on " `uname -a` #Only to display the machine where the job is running
echo "System release " `cat /etc/redhat-release` #And the system release
echo "CMSSW on Condor"

CMSSWVER=CMSSW_8_0_11
OUTDIR=root://cmseos.fnal.gov//store/user/rasharma/LHE_GEN_Analyzer_Output/
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_LO_NPle1_mjj200_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_LO_NPle1_DefaultCut_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_SM_mjj200_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_SM_DefaultCut_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_LO_NPle1_RunCardChanged_InitilaizeAllaQGCPar_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_ENuQQJJ_EWK_LO_SM_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_LO_NPle1_FT010e12_Pythia8CUEP8M1_13TeV_Madgraph.root
#OutPutFileName=aQGC_WPlepWMhadJJ_EWK_LO_NPle1_FT08e12_Pythia8CUEP8M1_13TeV_Madgraph.root
OutPutFileName=aQGC_WPlepWMhadJJ_EWK_LO_NPle1_FT014e12_Pythia8CUEP8M1_13TeV_Madgraph.root


echo ""
echo "parameter set:"
echo "CMSSWVER:   $CMSSWVER"

tar -xzf ${CMSSWVER}.tar.gz
cd ${CMSSWVER}
scram b ProjectRename
source /cvmfs/cms.cern.ch/cmsset_default.sh
# cmsenv
cd src/GEN-SIM-analyzer/GenAnalyzer
eval `scramv1 runtime -sh`

# run CMSSW
cmsRun python/ConfFile_cfg.py

# copy output to eos
echo "xrdcp output for condor"
for FILE in *.root
do
  #echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  #xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  echo "xrdcp -f ${FILE} ${OUTDIR}/${OutPutFileName}"
  xrdcp -f ${FILE} ${OUTDIR}/${OutPutFileName} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

echo ""
END_TIME=`date`
echo "finished at $END_TIME"
