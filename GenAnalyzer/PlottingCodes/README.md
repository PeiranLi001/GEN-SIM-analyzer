# Plot for aQGC different pars:
To get the aQGC distribution for different aQGC parameters on same canvas with ratio plot w.r.t. SM. 

	python Plot_aQGC.py


# List of available codes & uses

1. DataMCInfo.ini : This contains list of cross-section and number of events at MINIAOD level.

1. getBranchList.py : To print the list of all available branches in a tree

1. plot_functions.py : This contains many functions for plotting.

1. make_CompHistoFromSameRootFile.py : It uses plot_functions.py to compare two histogram from same file.

1. CompareBranchesFromSameTree.sh : Script to run make_CompHistoFromSameRootFile.py many times

1. compare_SM_with_RwgtSM_DiffPars.py : This macro compares SM with Reweighted SM for various initial parameters in customized cards.

1. Comp_Old_New_Sample.py : This macro compares inclusive polarized SM sample with the newly generated SM sample at generator level.

1. Signal_Reco_mT_plot.py : To plot a varialbe from MINIAOD signal sample after full event selection.

# To calculate cross-section for each parameter of aQGC

1. Get_Cross_Sec_aQGC.C  && Get_Cross_Sec_aQGC.h
2. Change input file name and output file name then run it.

# To plot aQGC par vs Cross-section

1. Make output of above step, input of macro ReadTextFile.py. This will create 9 different text file for each aQGC par

2. Then modify Template.py as per its structure

3. Create a cpp file using command:

		python makeCppGraphMaker.py

4. This will create a <File>.C. So, just run it like :

		root -l -b -q <File>.C
