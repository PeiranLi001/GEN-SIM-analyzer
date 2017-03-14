import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Options and Output Report
process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )



process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	#'file:/afs/cern.ch/user/r/rasharma/work/aQGC_Studies/MC_SampleGeneration/analyzeLHE/CMSSW_8_0_11/src/GenAnalyzer_Arun/genAnalyzer/test/SMP-RunIISummer15wmLHEGS-00029_142.root'
	'file:SMP-RunIISummer15wmLHEGS-00029_Hadronize_RunCard_ssWW.root'
	)
)

process.demo = cms.EDAnalyzer('GenAnalyzer',
	Verbose		=	cms.bool(False),

	genParticlesInputTag  = cms.InputTag("genParticles"),
	LHEEventInputTag = cms.InputTag("externalLHEProducer"),
	genJetsAK4jetsInputTag= cms.InputTag("ak4GenJets"),
	genJetsAK8jetsInputTag= cms.InputTag("ak8GenJets"),	#Not in use
	genMetTrueInputTag= cms.InputTag("genMetTrue"),		#not in use
	genMetCaloInputTag= cms.InputTag("genMetCalo")		#not in use

)

#process.GenAnalyzer = cms.EDProducer("GenAnalyzer",
#)

process.p = cms.Path(process.demo)
