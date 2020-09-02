# How To Use

```bash
cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src
cmsenv
git clone git@github.com:ram1123/GEN-SIM-analyzer.git
scramv1 b -j 8
cd GEN-SIM-analyzer
git submodule update --init
cd GenAnalyzer
# Befor running, check the config file and comment/uncomment line for MINIAOD/LHE-GEN
cmsRun python/ConfFile_cfg.py
```

# Submit job with [condor job](GenAnalyzer/condor_jobs)

1. To set proxy and create tarfile for your CMSSW directory:

   ```bash
	python create_cmssw_tarfile.py  
   ```

2. Then modify output root file name (or output directory) in script [RunGENSIM_condor.sh](GenAnalyzer/condor_jobs/RunGENSIM_condor.sh) and output log file name in script [RunGENSIM_condor.jdl](GenAnalyzer/condor_jobs/RunGENSIM_condor.jdl)

3. Then submit condor job using command:

   ```bash
	condor_submit RunGENSIM_condor.jdl
   ```

	
# Some General Script

1. Move crab files from crab directory hierchy to top directory : [Move_FileToTopDirFromCrab.py](GenAnalyzer/Move_FileToTopDirFromCrab.py)
2. To run cmsRun over more then 255 files: [GetListOfFile_More250files.py](GenAnalyzer/GetListOfFile_More250files.py)

# Plotting Macro

For instructions check the README inside the [plotting macro](GenAnalyzer/Plotting-Macro) directory.


# Plotting Macros (OLD):

You can find plotting macro [here](https://github.com/ram1123/GEN-SIM-analyzer/tree/LHEonlyMiniAOD/GenAnalyzer/PlottingCodes).

1. See the list of nevents, cross-section for all MiniAOD sampels: [DataMCInfo.ini](https://github.com/ram1123/GEN-SIM-analyzer/blob/LHEonlyMiniAOD/GenAnalyzer/PlottingCodes/DataMCInfo.ini)
2. To print the list of available branches: run script [getBranchList.py](GenAnalyzer/PlottingCodes_Old/getBranchList.py)
