#OutPutCodeName="aQGC_Plots_11Feb_Range_ModelFSonly_FS0_500em12_WPlepWMhad_Sensitivity"
#OutPutCodeName="aQGC_Plots_Range_FT0_12p5em12_WPlepWMhad_Sensitivity"
OutPutCodeName="aQGC_Plots_Range_FT0_12p5em12_WPlepWMhad"

datapath="/uscms_data/d3/rasharma/aQGC_analysis/AnalysisFramework/GENAnalyzer/LHEonlyGEN/CMSSW_8_0_11/src/GEN-SIM-analyzer/GenAnalyzer/PlottingCodes/CrossSec_Vs_aQGCpar/"

InputData=[ "FT0_12p5_fs0.txt", "FT0_12p5_fs1.txt", "FT0_12p5_ft0.txt", "FT0_12p5_ft1.txt", "FT0_12p5_ft2.txt", "FT0_12p5_fm0.txt", "FT0_12p5_fm1.txt", "FT0_12p5_fm6.txt", "FT0_12p5_fm7.txt"]
#InputData=[ "FT0_12p5_fs0_sensitivity.txt", "FT0_12p5_fs1_sensitivity.txt", "FT0_12p5_ft0_sensitivity.txt", "FT0_12p5_ft1_sensitivity.txt", "FT0_12p5_ft2_sensitivity.txt", "FT0_12p5_fm0_sensitivity.txt", "FT0_12p5_fm1_sensitivity.txt", "FT0_12p5_fm6_sensitivity.txt", "FT0_12p5_fm7_sensitivity.txt"]
	   

xlabel="aQGC Parameter (#times 10^{-12} GeV^{-4})"
ylabel="Cross-Section (pb)"
#ylabel="Sensitivity"


VarInTextFile=["aQGC_par","CrossSec","ErrCrossSec"]  # Each element of list corresponds to one column of intput text file

VarToPlot=[]

#GraphTitle="aQGC parameter Vs Sensitivity"
GraphTitle="aQGC parameter Vs Cross-Section"

tlatexx = [
	]

legends = [
	"FS0",	"FS1",
	"FT0",  "FT1",
	"FT2",
	"FM0",	"FM1",
	"FM6",	"FM7"
	]
	
iffit = 0
fitfunction = "pol3"
fitXrange = [-60.0,60.0]	# Comment this line if want to set auto range

yrange=[0.0,0.0]	# Y range; For auto range set yrange[0] = yrange[1]

xscaleFactor=0.0	# if you don't want to scale x-values then put it = 0.0
yscaleFactor=0.0	# if you don't want to scale y-values then put it = 0.0

pos1=0.15		# Legend x pos
pos2=0.83		# Legend y pos
sep=6.5			# y-distance b/w two legends fraction
cpos1=00.11		# CMS Prem x pos
cpos2=0.91		# CMS Prem y pos
cpos3=0.85		# DetInfo position
DetInfo="FT1"	# Detector Info


xoffset=1.2		# X-offset
yoffset=1.0		# Y-offset

setMaxdigit=4		# SetMaximumDigit
getstat=0		# GetStats
