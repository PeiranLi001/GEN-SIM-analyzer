# How To Use

	cmsrel CMSSW_8_0_11
	cd CMSSW_8_0_11/src
	cmsenv
	git clone git@github.com:ram1123/GEN-SIM-analyzer.git
	scramv1 b -j 8
	cd AnalyzeGenSim/GenAnalyzer
	cmsRun python/ConfFile_cfg.py
