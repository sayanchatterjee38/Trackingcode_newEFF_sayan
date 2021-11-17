if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    
    ### this is an example for running on RECO
    ### the name must be changed crab.cfg for actual running
    from CRABClient.UserUtilities import config
    config = config()

    config.General.requestName ='PbPb2018_thnsparse4dsumw2_phibinwidth0p1_Nov14_2021_Negtrk'
    config.General.workArea ='PbPb2018_thnsparse4dsumw2_phibinwidth0p1_Nov14_2021_Negtrk'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.psetName = 'run_PbPb_newcfg.py'
    config.JobType.pluginName = 'Analysis'
    
    
    config.Data.inputDBS = 'phys03'
    config.Data.inputDataset = '/MinBias_Hydjet_Drum5F_2018_5p02TeV/clindsey-RECODEBUG_20190625-5db5dfa073297cb96791f14c622e83e2/USER'
    config.Data.splitting = 'LumiBased' #'FileBased'
    config.Data.ignoreLocality = False
    config.Data.unitsPerJob = 15
    config.Data.totalUnits = -1
    config.Data.outLFNDirBase = '/store/user/sayan/' 
    config.Data.publication = False
    config.Data.outputDatasetTag = 'PbPb2018_thnsparse4dsumw2_phibinwidth0p1_Nov14_2021_Negtrk'
    
    #config.Site.storageSite = 'T2_IN_TIFR'
    config.Site.storageSite = 'T3_CH_CERNBOX'
    

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)
        
    
    submit(config)
