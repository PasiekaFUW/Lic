import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess


process = cms.Process('OmtfTree')

#
# For processing single files insert lines with 'file:/PATH/FILE.root'
# alernatively one can use 'root://xrootd.unl.edu//store/express/Run2015A.....root'
# or                       '/store/express/Run2015A/ExpressPhysics/FEVT/...root'
# (there is 255 file limit though). Can be empty for crab.
#


#dataDir='/scratch_cmsse/akalinow/CMS/Data/SingleMu/13_1_0_04_01_2024/SingleMu_ch0_OneOverPt_Run2029_13_1_0_04_01_2024/13_1_0_04_01_2024/240104_104036/0000/'
#GJ
dataDir='/scratch_cmsse/akalinow/CMS/Data/SingleMu/13_1_0_11_03_2024/SingleMu_ch0_OneOverPt_Run2023_13_1_0_11_03_2024/13_1_0_11_03_2024/240311_101428/0000/'
lsCommand='ls -1 '+dataDir+'|grep root | grep m_99'

print ('command: ',lsCommand)
lsOutput= subprocess.Popen(lsCommand, stdout=subprocess.PIPE, shell=True, text=True).communicate()[0]
files=[]
for f in lsOutput.split():
  files.append('file:'+dataDir+f)
print ("Number of files in direcotry:",dataDir," ---> ", len(files))


process.source = cms.Source("PoolSource", 
#GJ fileNames = cms.untracked.vstring( 'file:data.root',),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/OMTF/PrivateProductionForOMTFStudy/13_1_0_03_04_2024/SingleMu_ch0_OneOverPt_Run2029_13_1_0_03_04_2024/13_1_0_03_04_2024/240403_080928/0000/SingleMu_OneOverPt_1_100_m_1.root'),      
fileNames = cms.untracked.vstring(files),
#skipEvents =  cms.untracked.uint32(220)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(44000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(44) ) #debug

#
# import of standard configurations
#
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryDB_cff')
#process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
#process.load('Configuration.Geometry.GeometryExtendedRun4D95Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


#
# set proper GlobalTag
#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '130X_dataRun3_Prompt_v3'
#process.GlobalTag.globaltag = '141X_mcRun4_realistic_v3'
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')


#
# message logger
#
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.enableStatistics = False
process.MessageLogger.cout.enable =False 
process.MessageLogger.suppressWarning  = cms.untracked.vstring('*')

# Calibrate Digis
process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
process.CalibratedDigis.dtDigiTag = "simMuonDTDigis" 
process.CalibratedDigis.scenario = 0

# DTTriggerPhase2
process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.dtTriggerPhase2PrimitiveDigis.debug = False
process.dtTriggerPhase2PrimitiveDigis.dump = False
process.dtTriggerPhase2PrimitiveDigis.scenario = 0
#Warto dodac ponizsze co wylacza "coincidence filter". Jak to dizala - MK?
#process.dtTriggerPhase2PrimitiveDigis.co_option = -1

process.licDigiAnalysis = cms.EDAnalyzer("LicDigiAnalysis",

  srcDTPh_leg = cms.InputTag('simDtTriggerPrimitiveDigis'),
  srcDTTh_leg = cms.InputTag('simDtTriggerPrimitiveDigis'),
  srcDTPh_upg = cms.InputTag('dtTriggerPhase2PrimitiveDigis'),
  srcDTTh_upg = cms.InputTag('dtTriggerPhase2PrimitiveDigis'),
)

process.out = cms.OutputModule("PoolOutputModule", 
   fileName = cms.untracked.string("superprimitives2.root"),
   outputCommands = cms.untracked.vstring(
                        "drop *",
                        "keep *_*_*_OmtfTree",

     )  
)

process.digi_step =  cms.Path(process.CalibratedDigis * process.dtTriggerPhase2PrimitiveDigis * process.licDigiAnalysis)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.digi_step, process.endjob_step)


#process.output_step=cms.EndPath(process.out)
#process.schedule.extend([process.output_step])

#print process.dumpPython();
#print process.schedule
