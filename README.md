# How To Use

	cmsrel CMSSW_8_0_11
	cd CMSSW_8_0_11/src
	cmsenv
	git clone git@github.com:ram1123/GEN-SIM-analyzer.git
	scramv1 b -j 8
	cd AnalyzeGenSim/GenAnalyzer
	# Befor running, check the config file and comment/uncomment line for MINIAOD/LHE-GEN
	cmsRun python/ConfFile_cfg.py

* Also, you need to add the input root file. If you have long list of root file then you can use the script [GetListOfFile_More250files.py](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/GetListOfFile_More250files.py)
* Here in the code you need to provide the path of input root files.
* Copy the output of the above code and paste it inot python/ConfFile_cfg.py

## Submit job with [condor job](https://github.com/ram1123/GEN-SIM-analyzer/tree/LHEonlyMiniAOD/GenAnalyzer/condor_jobs)

1. To set proxy and create tarfile for your CMSSW directory:
	
		source create_cmssw_tarfile.sh  # bash shell

2. Then modify output root file name (or output directory) in script [RunGENSIM_condor.sh](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/condor_jobs/RunGENSIM_condor.sh) and output log file name in script [RunGENSIM_condor.jdl](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/condor_jobs/RunGENSIM_condor.jdl)

3. Then submit condor job using command:

		condor_submit RunGENSIM_condor.jdl

	


# Some General Script

1. Move crab files from crab directory hierchy to top directory : [Move_FileToTopDirFromCrab.py](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/Move_FileToTopDirFromCrab.py)
2. To run cmsRun over more then 255 files: [GetListOfFile_More250files.py](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/GetListOfFile_More250files.py)



# Plotting Macros:

You can find plotting macro [here](https://github.com/ram1123/GEN-SIM-analyzer/tree/LHEonlyMiniAOD/GenAnalyzer/PlottingCodes).

1. See the list of nevents, cross-section for all MiniAOD sampels: [DataMCInfo.ini](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/PlottingCodes/DataMCInfo.ini)
2. To print the list of available branches: run script [getBranchList.py](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/PlottingCodes/getBranchList.py)
