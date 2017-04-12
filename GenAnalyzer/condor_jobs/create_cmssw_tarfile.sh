#eval export PATH="/cvmfs/cms.cern.ch/share/overrides/bin:/uscms_data/d3/rasharma/aQGC_analysis/AnalysisFramework/GENAnalyzer/LHEonlyGEN/CMSSW_8_0_11/bin/slc6_amd64_gcc530:/uscms_data/d3/rasharma/aQGC_analysis/AnalysisFramework/GENAnalyzer/LHEonlyGEN/CMSSW_8_0_11/external/slc6_amd64_gcc530/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_11/bin/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_11/external/slc6_amd64_gcc530/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.7.1-giojec/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin:/usr/krb5/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin";
eval `scramv1 runtime -sh`
PCHECK=`voms-proxy-info -timeleft`
if [[ "$PCHECK" -eq 0 ]]; then
  voms-proxy-init -voms cms --valid 168:00
fi
tar  -zcf CMSSW_8_0_11.tar.gz -C /uscms_data/d3/rasharma/aQGC_analysis/AnalysisFramework/GENAnalyzer/LHEonlyGEN/CMSSW_8_0_11/.. CMSSW_8_0_11
