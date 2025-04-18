import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess

#PPG
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

#process = cms.Process('OmtfTree')
#PPG
process = cms.Process('Analysis',Phase2C17I13M9)

#
# For processing single files insert lines with 'file:/PATH/FILE.root'
# alernatively one can use 'root://xrootd.unl.edu//store/express/Run2015A.....root'
# or                       '/store/express/Run2015A/ExpressPhysics/FEVT/...root'
# (there is 255 file limit though). Can be empty for crab.
#

#Sample used up to zad 7
#dataDir='/scratch_cmsse/akalinow/CMS/Data/SingleMu/13_1_0_04_01_2024/SingleMu_ch0_OneOverPt_Run2029_13_1_0_04_01_2024/13_1_0_04_01_2024/240104_104036/0000/'
#Sample used from zad 7 including
dataDir='/scratch_cmsse/akalinow/CMS/Data/SingleMu/13_1_0_11_03_2024/SingleMu_ch0_OneOverPt_Run2023_13_1_0_11_03_2024/13_1_0_11_03_2024/240311_101428/0000/'
#Sample used for propagation not working due to lacking objects of PSimHit
#dataDir='/scratch_cmsse/akalinow/CMS/Data/SingleMu/13_1_0_13_02_2024/SingleMu_ch0_OneOverPt_Run2029_13_1_0_13_02_2024/13_1_0_13_02_2024/240213_094022/0000/'
#PPG
#dataDir='/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/OMTF/13_1_0_11_03_2024/SingleMu_ch0_OneOverPt_Run2023_13_1_0_11_03_2024/13_1_0_11_03_2024/240311_101428/0000/'
lsCommand='ls -1 '+dataDir+'|grep root | grep m_99' 

print ('command: ',lsCommand)
lsOutput= subprocess.Popen(lsCommand, stdout=subprocess.PIPE, shell=True, text=True).communicate()[0]
files=[]
for f in lsOutput.split():
  files.append('file:'+dataDir+f)
print ("Number of files in direcotry:",dataDir," ---> ", len(files))

process.source = cms.Source("PoolSource", 
#fileNames = cms.untracked.vstring(
# 'root://xrootd-cms.infn.it//store/data/Run2023B/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/367/079/00000/b01b794e-9075-4d42-b273-85b2bc66f13a.root',
#  ),
fileNames = cms.untracked.vstring(files),
#skipEvents =  cms.untracked.uint32(220)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(44000) ) #Number of Events, Max == 44000

#
# import of standard configurations
''' Before #PPG
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
'''
#PPG
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D99Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')





#
# set proper GlobalTag
#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '130X_dataRun3_Prompt_v3'
#PPG
process.GlobalTag.globaltag = '131X_mcRun4_realistic_v7'

#
# message logger
#
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 4000 #mozna modyfikowac 
process.MessageLogger.cerr.enableStatistics = False
process.MessageLogger.cout.enable =False 
process.MessageLogger.suppressWarning  = cms.untracked.vstring('*')

process.licDigiAnalysis = cms.EDAnalyzer("LicDigiAnalysis",

  srcDTPh_leg = cms.InputTag('simDtTriggerPrimitiveDigis'),
  srcDTTh_leg = cms.InputTag('simDtTriggerPrimitiveDigis'),
)


process.digi_step =  cms.Path(process.licDigiAnalysis)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.digi_step, process.endjob_step)


#print process.dumpPython();
#print process.schedule
