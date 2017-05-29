import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Options and Output Report
process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	#'file:/afs/cern.ch/user/r/rasharma/work/aQGC_Studies/MC_SampleGeneration/analyzeLHE/CMSSW_8_0_11/src/GenAnalyzer_Arun/genAnalyzer/test/SMP-RunIISummer15wmLHEGS-00029_142.root'
	#'file:/uscms_data/d3/rasharma/aQGC_analysis/AnalysisFramework/GENAnalyzer/CMSSW_8_0_11/src/GEN-SIM-analyzer/GenAnalyzer/MiniAOD_00E2D4C8-10C8-E611-AF91-B8CA3A70A5E8.root'
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15wmLHEGS/WWJJ_SS_WToLNu_EWK_aQGC-FT-FS-FM_TuneCUETP8M1_13TeV_madgraph-pythia8/LHE/MCRUN2_71_V1-v1/80000/60BA27D6-2BEC-E611-AC06-003048CDCDD2.root'
	'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWJJ_SS_WToLNu_EWK_aQGC-FT-FS-FM_TuneCUETP8M1_13TeV_madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/58D945A3-60F1-E611-9F9C-00259048BF92.root'
		),
    skipBadFiles = cms.untracked.bool(True)
)


process.demo = cms.EDAnalyzer('GenAnalyzer',
	Verbose		=	cms.bool(False),
	#genParticlesInputTag  = cms.InputTag("prunedGenParticles"),	# Uncomment if running on MiniAOD
	#LHEEventInputTag = cms.InputTag("source"),			# Uncomment if running on MiniAOD
	genParticlesInputTag  = cms.InputTag("genParticles"),		# Uncomment if running on GEN only sample
	LHEEventInputTag = cms.InputTag("externalLHEProducer"),		# Uncomment if running on GEN only

)

#process.GenAnalyzer = cms.EDProducer("GenAnalyzer",
#)

process.p = cms.Path(process.demo)
