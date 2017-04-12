import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Options and Output Report
process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	#'file:/afs/cern.ch/user/r/rasharma/work/aQGC_Studies/MC_SampleGeneration/analyzeLHE/CMSSW_8_0_11/src/GenAnalyzer_Arun/genAnalyzer/test/SMP-RunIISummer15wmLHEGS-00029_142.root'
	#'file:/uscms_data/d3/rasharma/aQGC_analysis/AnalysisFramework/GENAnalyzer/CMSSW_8_0_11/src/GEN-SIM-analyzer/GenAnalyzer/MiniAOD_00E2D4C8-10C8-E611-AF91-B8CA3A70A5E8.root'
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWJJToLNuQQ_LL_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/128AFDAD-C7C6-E611-A609-001E67A3F49D.root'
       ' root://cmseos.fnal.gov//store/user/rasharma/CMSSW_FullSimulation/aQGC_WPlepWMhadJJ_EWK_DefaultCut_Pythia8CUEP8M1_13TeV_Madgraph/RunIISummer15wmLHEGS/170403_222535/0000/SMP-RunIISummer15wmLHEGS-00029_Hadronize_MadDefCard_106.root',
	),
    skipBadFiles = cms.untracked.bool(True)
)


process.demo = cms.EDAnalyzer('GenAnalyzer',
	Verbose		=	cms.bool(False),

	LHEEventInputTag = cms.InputTag("source"),	# Uncomment this if running on MiniAOD
	#LHEEventInputTag = cms.InputTag("externalLHEProducer"),	# Uncomment this if running on GEN only

)

#process.GenAnalyzer = cms.EDProducer("GenAnalyzer",
#)

process.p = cms.Path(process.demo)
