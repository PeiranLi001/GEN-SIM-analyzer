import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# Options and Output Report
process.options = cms.untracked.PSet(
                                     wantSummary = cms.untracked.bool(True)
                                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring(FileUtils.loadListFromFile('GF_HH_Benchmark2.txt'))
#readFiles = cms.untracked.vstring(FileUtils.loadListFromFile('TEMP_NAME.txt'))


process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            #fileNames = cms.untracked.vstring(
                            ##'file:/afs/cern.ch/user/r/rasharma/work/aQGC_Studies/MC_SampleGeneration/analyzeLHE/CMSSW_8_0_11/src/GenAnalyzer_Arun/genAnalyzer/test/SMP-RunIISummer15wmLHEGS-00029_142.root'
                            #	),
                            fileNames = readFiles,
                            skipBadFiles = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )

process.demo = cms.EDAnalyzer('GenAnalyzer',
                              Verbose		=	cms.bool(False),
                              OutPutFileName = cms.string("test.root"),
                              #OutPutFileName = cms.string("GF_HH_Benchmark222.root"),
                              #genParticlesInputTag  = cms.InputTag("prunedGenParticles"),	# Uncomment if running on MiniAOD
                              #LHEEventInputTag = cms.InputTag("source"),			# Uncomment if running on MiniAOD
                              genParticlesInputTag  = cms.InputTag("genParticles"),		# Uncomment if running on GEN only sample
                              LHEEventInputTag = cms.InputTag("externalLHEProducer"),		# Uncomment if running on GEN only
                              ak4GenJetsInputTag = cms.InputTag("ak4GenJets"),
                              ak8GenJetsInputTag = cms.InputTag("ak8GenJets"),
                              
                              )

#process.GenAnalyzer = cms.EDProducer("GenAnalyzer",
#)

process.p = cms.Path(process.demo)
