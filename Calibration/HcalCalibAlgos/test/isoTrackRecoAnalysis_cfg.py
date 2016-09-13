import FWCore.ParameterSet.Config as cms
process = cms.Process('ANALYSIS')

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag=autoCond['run2_mc']

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('warnings','errors',
                                         'cout','cerr'),
    categories = cms.untracked.vstring('HcalIsoTrack'), 
    debugModules = cms.untracked.vstring('*'),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    cerr = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        noTimeStamps = cms.untracked.bool(False),
        FwkReport = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(10000),
            limit = cms.untracked.int32(10000000)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        FwkJob = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        FwkSummary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(10000),
            limit = cms.untracked.int32(10000000)
        ),
        threshold = cms.untracked.string('INFO')
     ),
     cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        noTimeStamps = cms.untracked.bool(True),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        HcalIsoTrack = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
       )
    )
)

process.load('Calibration.HcalCalibAlgos.isoAnalyzer_m2_cfi')
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
#        'file:root://eoscms//eos/cms//store/user/qnguyen/HI2013_Run211313/crab_HI2013/160826_121802/0000/DataRereco_1.root',
#        'file:root://eoscms//eos/cms//store/user/qnguyen/HI2013_Run211313/crab_HI2013/160826_121802/0000/DataRereco_2.root',
#        'file:root://eoscms//eos/cms//store/user/qnguyen/HI2013_Run211313/crab_HI2013/160826_121802/0000/DataRereco_3.root',
#        'file:root://eoscms//eos/cms//store/user/qnguyen/HI2013_Run211313/crab_HI2013/160826_121802/0000/DataRereco_4.root',
#        'file:root://eoscms//eos/cms//store/user/qnguyen/HI2013_Run211313/crab_HI2013/160826_121802/0000/DataRereco_5.root'
'root://eoscms.cern.ch//store/group/dpg_hcal/comm_hcal/RecoAlgos/HIN2013/MinimumBias/160830_155324/DataRereco_1.root',
'root://eoscms.cern.ch//store/group/dpg_hcal/comm_hcal/RecoAlgos/HIN2013/MinimumBias/160830_155324/DataRereco_10.root',
'root://eoscms.cern.ch//store/group/dpg_hcal/comm_hcal/RecoAlgos/HIN2013/MinimumBias/160830_155324/DataRereco_101.root',
'root://eoscms.cern.ch//store/group/dpg_hcal/comm_hcal/RecoAlgos/HIN2013/MinimumBias/160830_155324/DataRereco_102.root',
'root://eoscms.cern.ch//store/group/dpg_hcal/comm_hcal/RecoAlgos/HIN2013/MinimumBias/160830_155324/DataRereco_103.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_1.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_2.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_3.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_4.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_5.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_6.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_7.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_8.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_9.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_X_2016-09-07-1100/src/Phase1/step3_phase1_10.root',
#'file:root://eoscms.cern.ch//store/user/qnguyen/MinimumBias/crab_HI2013_Json/160830_155324/0000/DataRereco_139.root',
#'file:/afs/cern.ch/work/q/qnguyen/public/PulseShape/CMSSW_8_1_0_pre10/src/MCGenerator/Phase1/step3_phase1.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.TFileService = cms.Service("TFileService",
   fileName = cms.string('IsoTrack.root')
)

process.p = cms.Path(process.HcalIsoTrkAnalyzer)

