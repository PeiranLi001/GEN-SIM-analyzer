import os
import ROOT as ROOT
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_SMNoDecay_mjj200_Pythia8CUEP8M1_13TeV_Madgraph/RunIISummer15wmLHEGS/170410_120247/'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_mjj200_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_DefaultCut_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_SM_mjj200_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_SM_DefaultCut_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_RunCardChanged_InitilaizeAllaQGCPar_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_ENuQQJJ_EWK_LO_SM_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_FT010e12_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_FT08e12_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_FT014e12_Pythia8CUEP8M1_13TeV_Madgraph'
#source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_SM_Pythia8CUEP8M1_13TeV_Madgraph_RwgtValidation_v2'
source = '/eos/uscms/store/user/rasharma/CMSSW_FullSimulation_April2017/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_Pythia8CUEP8M1_13TeV_Madgraph_RwgtValidation'

search1='log'
search2='failed'
Arrayfilepath = []
#for root, dirs, filenames in os.walk(source):
for root, dirs, filenames in os.walk(source):	
	for f in filenames:
		filepath = root + os.sep + f
		#print f
		if filepath.find(search2) == -1:	# Don't select file if it contains word failed in path
			if filepath.find(search1) == -1:	# Don't select file if it contains wored log in path
				if filepath.endswith(".root"):	# select only root files
					Arrayfilepath.append(filepath)
#print len(Arrayfilepath)					
#print Arrayfilepath[0]

"""
for i in range(0,len(Arrayfilepath)):
	file = ROOT.TFile(Arrayfilepath[i])
	print Arrayfilepath[i]
	if not file:
		print 'Failed to open %s' % Arrayfilepath[i]
		exit(0)
	tree = file.Get('Events')
	print tree.GetEntries()
"""	
print "myfilelist = cms.untracked.vstring()"
print "myfilelist.extend( ["
for i in range(0,len(Arrayfilepath)):
	NewFilePath = Arrayfilepath[i].replace("/eos/uscms","")
	if ((i+1)%250==0):
		print "\t'root://cmseos.fnal.gov/"+NewFilePath+"'"
		print "])"
		print "myfilelist.extend( ["
	else:
		print "\t'root://cmseos.fnal.gov/"+NewFilePath+"',"
print "])"
print "process.source = cms.Source(\"PoolSource\","
print "\t\tfileNames = myfilelist,"
print "\t\tskipBadFiles = cms.untracked.bool(True)"
print "\t\t)"
