from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'Signal_LL'
config.General.workArea = 'Signal_LL'
config.section_('JobType')
config.JobType.psetName = 'python/ConfFile_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.pyCfgParams = ['global_tag=80X_mcRun2_asymptotic_2016_TrancheIV_v7','leptonFilter=False', 'MC=True', 'isCrab=True', 'DoJECCorrection=True', 'DoPuppi=True','ReDoPruningAndSoftdropPuppi=True']
#config.JobType.inputFiles = ['Spring16_23Sep2016V2_MC_L1FastJet_AK8PFchs.txt','Spring16_23Sep2016V2_MC_L2Relative_AK8PFchs.txt','Spring16_23Sep2016V2_MC_L3Absolute_AK8PFchs.txt','Spring16_23Sep2016V2_MC_L1FastJet_AK4PFchs.txt','Spring16_23Sep2016V2_MC_L2Relative_AK4PFchs.txt','Spring16_23Sep2016V2_MC_L3Absolute_AK4PFchs.txt','Spring16_23Sep2016V2_MC_Uncertainty_AK4PFchs.txt','Spring16_23Sep2016V2_MC_Uncertainty_AK8PFchs.txt','Spring16_23Sep2016V2_MC_L1FastJet_AK8PFPuppi.txt','Spring16_23Sep2016V2_MC_L2Relative_AK8PFPuppi.txt','Spring16_23Sep2016V2_MC_L3Absolute_AK8PFPuppi.txt','Spring16_23Sep2016V2_MC_L1FastJet_AK4PFPuppi.txt','Spring16_23Sep2016V2_MC_L2Relative_AK4PFPuppi.txt','Spring16_23Sep2016V2_MC_L3Absolute_AK4PFPuppi.txt','Spring16_23Sep2016V2_MC_Uncertainty_AK4PFPuppi.txt','Spring16_23Sep2016V2_MC_Uncertainty_AK8PFPuppi.txt' ]
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['LHEinfo.root']
#config.JobType.maxMemoryMB = 2500    # 2.5 GB     
config.JobType.maxJobRuntimeMin = 900 #15 h
config.section_('Data')
config.Data.inputDataset = '/WWJJToLNuQQ_LL_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.unitsPerJob = 1
config.Data.inputDBS = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/rasharma/LHE_GEN_Analyzer_Output_13April/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.blacklist= ['T2_US_Purdue','T2_UA_KIPT']

#NB: SAMPLES HAVE TO BE UPDATED!

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand

    #Make sure you set this parameter (here or above in the config it does not matter)
    config.General.workArea = 'Crab_Logs_25April2017'

    def submit(config):
        res = crabCommand('submit', config = config)

    #########    From now on that's what users should modify: this is the a-la-CRAB2 configuration part.

   
    config.General.requestName = 'ssWW_aQGC_EWK'
    config.Data.inputDataset = '/WWJJ_SS_WToLNu_EWK_aQGC-FT-FS-FM_TuneCUETP8M1_13TeV_madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.Data.outLFNDirBase = '/store/user/rasharma/LHE_GEN_Analyzer_Output_13April/ssWW_aQGC_EWK/'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'ssWW_SM_EWK'
    config.Data.inputDataset = '/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.Data.outLFNDirBase = '/store/user/rasharma/LHE_GEN_Analyzer_Output_13April/ssWW_SM_EWK/'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
