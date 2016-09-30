import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                                  '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/08080C9F-C078-E611-AEF1-FA163E5647FC.root',
                                  '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/084B4FBE-C078-E611-9740-02163E013993.root',
                                  '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/0A6B1BD0-C078-E611-A332-02163E014113.root'
                        ),
                    secondaryFileNames = cms.untracked.vstring(),
#                   lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v11'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


basePath = '/afs/cern.ch/user/d/deguio/public/qualityPerSub_2016/'
subsystemList = ['Pix','Strip','Ecal','Hcal','Dt','Rpc','Es','Csc','Track']
fileList = []
for sub in subsystemList:
    fileList.append(basePath+'Cert_13TeV_2016_'+sub+'.txt')


process.MyAnalysis =cms.EDAnalyzer("MyMiniAODAnalyzer",

                       #caloJetTag              = cms.untracked.InputTag("slimmedJets"),
                       PFJetTag                = cms.untracked.InputTag("slimmedJets"),
                       vtx                     = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                                   
                       maxJetEta               = cms.untracked.double(3.0),
                       minJetPt                = cms.untracked.double(30.0),

                       lumiFile                = cms.untracked.string(basePath+'run_ls_lumi_2016.txt'),
                       subsystems              = cms.untracked.vstring(subsystemList),
                       qualityFiles            = cms.untracked.vstring(fileList),
                       quantiles               = cms.untracked.vdouble(0.,0.25,0.50,0.75,1.)

                       )


process.mypath  = cms.Path(process.MyAnalysis)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("miniAODTree.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )
 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options =  cms.untracked.PSet(
                   #allowUnscheduled = cms.untracked.bool(True),
                   wantSummary = cms.untracked.bool(True),
                   )


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
    debugModules  = cms.untracked.vstring('*')
)
