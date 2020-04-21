import sys
sys.path.append('..')
import json
import os
from itertools import combinations
import onPPILinkPred as ppiLPred
import bioGRID as bg
import IntAct
import STRING as string
import numpy as np
import traversalHelper as tr
import pandas as pd
import domain as dm
from collections import defaultdict
import time
import MINT

class helper:
    @staticmethod
    def write_runTime(tag, runTime):
        fName = './resultData/runTime.json'
        if not os.path.exists(fName):
            with open(fName, 'w') as f: f.write(json.dumps({}))
        
        runTimeD = {}
        with open(fName, 'r') as f: runTimeD = json.loads(f.read())
        runTimeD[tag] = runTime
        with open(fName, 'w') as f: f.write(json.dumps(runTimeD))

    @staticmethod
    def trim_ppi_result(fNames, datasetClass):
        # get only the top PPI & scores equal to the size of its original dataset
        trimNum = {'bioGRID': int(len([*bg.parse_bioGRID(root='../')][1].index)*0.5)
        , 'STRING': int(len([*string.parse_STRING(root='../')][1].index)*0.5)
        , 'IntAct': int(len([*IntAct.parse_IntAct(root='../')][1].index)*0.5)
        , 'MINT': int(len([*MINT.parse_MINT(root='../')][1].index)*0.5)
        , 'IntAct_spoke': int(len([*IntAct.parse_IntAct(root='../', spokeModel=True)][1].index)*0.5)}

        if not os.path.exists('./resultData/trimmed_predPPIs.json'):
            with open('./resultData/trimmed_predPPIs.json', 'w') as f: pass
            with open('./resultData/trimmed_predScores.json', 'w') as f: pass
        
        for i in range(len(fNames)):
            predPPI, predScore = [], []
            with open("./resultData/{}_PPI.json".format(fNames[i]), 'r') as f:
                for line in f.readlines(): predPPI.append(json.loads(line)[0:trimNum[datasetClass[i]]])
            with open("./resultData/{}_score.json".format(fNames[i]), 'r') as f:
                for line in f.readlines(): predScore.append(json.loads(line)[0:trimNum[datasetClass[i]]])
            predPPIs, predScores = {fNames[i]: predPPI}, {fNames[i]: predScore}
            with open('./resultData/trimmed_predPPIs.json', 'a+') as f: f.write(json.dumps(predPPIs)+"\n")
            with open('./resultData/trimmed_predScores.json', 'a+') as f: f.write(json.dumps(predScores)+"\n")

    @staticmethod
    def trim_multiple_ppi_result(fNames, datasetClass, trialSize):
        trimNum = {'bioGRID': int(len([*bg.parse_bioGRID(root='../')][1].index)*0.5)
        , 'STRING': int(len([*string.parse_STRING(root='../')][1].index)*0.5)
        , 'IntAct': int(len([*IntAct.parse_IntAct(root='../')][1].index)*0.5)
        , 'MINT': int(len([*MINT.parse_MINT(root='../')][1].index)*0.5)
        , 'IntAct_spoke': int(len([*IntAct.parse_IntAct(root='../', spokeModel=True)][1].index)*0.5)}

        if not os.path.exists('./resultData/trimmed_predPPIs.json'):
            with open('./resultData/trimmed_predPPIs.json', 'w') as f: pass
            with open('./resultData/trimmed_predScores.json', 'w') as f: pass
        
        for i in range(len(fNames)):
            predPPI, predScore = [], []
            for j in range(trialSize):
                with open("./resultData/{}_{}_PPI.json".format(fNames[i], j), 'r') as f:
                    predPPI.append(json.loads(f.read())[0:trimNum[datasetClass[i]]])
                with open("./resultData/{}_{}_score.json".format(fNames[i], j), 'r') as f:
                    predScore.append(json.loads(f.read())[0:trimNum[datasetClass[i]]])
            predPPIs, predScores = {fNames[i]: predPPI}, {fNames[i]: predScore}
            with open('./resultData/trimmed_predPPIs.json', 'a+') as f: f.write(json.dumps(predPPIs)+"\n")
            with open('./resultData/trimmed_predScores.json', 'a+') as f: f.write(json.dumps(predScores)+"\n")

    @staticmethod
    def write_resultData(predPPI, predScore, fName):
        with open('./resultData/{}_PPI.json'.format(fName), 'w') as f:
            f.write(json.dumps(predPPI))
        with open('./resultData/{}_score.json'.format(fName), 'w') as f:
            f.write(json.dumps(predScore))

    @staticmethod
    def append_precRecMap_multiCore(fNames, predPPI, samplePPI, datasetClass, coreNo, isGGI=False, logging=False):
        if isGGI: i = 0
        else: i = 1
        fullPPISet = {'bioGRID': [list(ppi) for ppi in np.asarray([*bg.parse_bioGRID(root='../')][i][['nodeA', 'nodeB']])]
        , 'STRING': [list(ppi) for ppi in np.asarray([*string.parse_STRING(root='../')][i][['nodeA', 'nodeB']])]
        , 'IntAct': [list(ppi) for ppi in np.asarray([*IntAct.parse_IntAct(root='../')][i][['nodeA', 'nodeB']])]
        , 'MINT': [list(ppi) for ppi in np.asarray([*MINT.parse_MINT(root='../')][i][['nodeA', 'nodeB']])]
        , 'IntAct_spoke': [list(ppi) for ppi in np.asarray([*IntAct.parse_IntAct(root='../', spokeModel=True)][i][['nodeA', 'nodeB']])]}

        fullPrecRecMap = {}
        if not os.path.exists("./resultData/PRCurveMap.json"):
            with open("./resultData/PRCurveMap.json", "w") as f: f.write(json.dumps(fullPrecRecMap))

        precRecMap = ppiLPred.precRecMap_multiCore(fNames, predPPI, samplePPI, [fullPPISet[i] for i in datasetClass], coreNo, logging)
        with open('./resultData/PRCurveMap.json', 'r') as f: fullPrecRecMap = json.loads(f.read())
        fullPrecRecMap.update(precRecMap)
        with open('./resultData/PRCurveMap.json', 'w') as f: f.write(json.dumps(fullPrecRecMap))


class initDataGen:
    @staticmethod
    def CN_tenTrial():
        tags = ['commonNeighbor']

        #bioGRID
        bioGRID_samplePPIs = []
        with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
        for t in range(len(tags)):
            for i in range(len(bioGRID_samplePPIs)):
                bioGRID_samplePPI = bioGRID_samplePPIs[i]
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(bioGRID_samplePPI, tags[t], scoreArgsDict={}, coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime('{}_tenTrial_bioGRID_{}'.format(tags[t], i), runTime)
                helper.write_resultData(predPPI, predScore, '{}_tenTrial_bioGRID_{}'.format(tags[t], i))

    @staticmethod
    def L3_uvJoin_tenTrial():
        #bioGRID
        bioGRID_samplePPIs = []
        with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
        for i in range(len(bioGRID_samplePPIs)):
            bioGRID_samplePPI = bioGRID_samplePPIs[i]
            startTime = time.time()
            predPPI, predScore = ppiLPred.multiCore_PPILinkPred(bioGRID_samplePPI, 'L3uvJoin', scoreArgsDict={}, coreNo=14, logging=True)
            runTime = time.time()-startTime
            helper.write_runTime('L3uvJoin_tenTrial_bioGRID_{}'.format(i), runTime)
            helper.write_resultData(predPPI, predScore, 'L3uvJoin_tenTrial_bioGRID_{}'.format(i))

    @staticmethod
    def ExactL3_uvJoin_tenTrial():
        tags = ['xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        scoreDicts = [{'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}, {'xyContrib': 'basic', 'uvSpec': 'basic', 'xySpec': 'basic', 'uvJoin': 'basic'}]

        #bioGRID
        bioGRID_samplePPIs = []
        with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
        for t in range(len(tags)):
            for i in range(len(bioGRID_samplePPIs)):
                fName = '{}_tenTrial_bioGRID_{}'.format(tags[t], i)
                bioGRID_samplePPI = bioGRID_samplePPIs[i]
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(bioGRID_samplePPI, 'interStr', scoreArgsDict=scoreDicts[t], coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def MINT_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        methods = ['commonNeighbor', 'L3uvJoin', 'interStr', 'interStr']
        scoreArgDicts = [{}, {}, {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}, {'xyContrib': 'basic', 'uvSpec': 'basic', 'xySpec': 'basic', 'uvJoin': 'basic'}]

        MINT_samplePPIs = []
        with open("./genData/MINT_sampledPPIs.json", "r") as f: MINT_samplePPIs = json.loads(f.read())

        for t in range(len(tags)):
            for i in range(len(MINT_samplePPIs)):
                fName = '{}_tenTrial_MINT_{}'.format(tags[t], i)
                MINT_samplePPI = MINT_samplePPIs[i]
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(MINT_samplePPI, methods[t], scoreArgsDict=scoreArgDicts[t], coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def IntAct_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        methods = ['commonNeighbor', 'L3uvJoin', 'interStr', 'interStr']
        scoreArgDicts = [{}, {}, {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}, {'xyContrib': 'basic', 'uvSpec': 'basic', 'xySpec': 'basic', 'uvJoin': 'basic'}]

        IntAct_samplePPIs = []
        with open("./genData/IntAct_sampledPPIs.json", "r") as f: IntAct_samplePPIs = json.loads(f.read())

        for t in range(len(tags)):
            for i in range(len(IntAct_samplePPIs)):
                fName = '{}_tenTrial_IntAct_{}'.format(tags[t], i)
                IntAct_samplePPI = IntAct_samplePPIs[i]
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(IntAct_samplePPI, methods[t], scoreArgsDict=scoreArgDicts[t], coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def IntAct_spoke_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        methods = ['commonNeighbor', 'L3uvJoin', 'interStr', 'interStr']
        scoreArgDicts = [{}, {}, {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}, {'xyContrib': 'basic', 'uvSpec': 'basic', 'xySpec': 'basic', 'uvJoin': 'basic'}]

        IntAct_samplePPIs = []
        with open("./genData/IntAct_spoke_sampledPPIs.json", "r") as f: IntAct_samplePPIs = json.loads(f.read())

        for t in range(len(tags)):
            for i in range(len(IntAct_samplePPIs)):
                fName = '{}_tenTrial_IntAct_spoke_{}'.format(tags[t], i)
                IntAct_samplePPI = IntAct_samplePPIs[i]
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(IntAct_samplePPI, methods[t], scoreArgsDict=scoreArgDicts[t], coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def STRING_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        methods = ['commonNeighbor', 'L3uvJoin', 'interStr', 'interStr']
        scoreArgDicts = [{}, {}, {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}, {'xyContrib': 'basic', 'uvSpec': 'basic', 'xySpec': 'basic', 'uvJoin': 'basic'}]

        STRING_samplePPIs = []
        with open("./genData/STRING_sampledPPIs.json", "r") as f: STRING_samplePPIs = json.loads(f.read())

        for t in range(len(tags)):
            for i in range(len(STRING_samplePPIs)):
                fName = '{}_tenTrial_STRING_{}'.format(tags[t], i)
                STRING_samplePPI = STRING_samplePPIs[i]
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(STRING_samplePPI, methods[t], scoreArgsDict=scoreArgDicts[t], coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def CRA_tenTrial():
        for dataset in ['bioGRID', 'STRING', 'MINT']:
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f:
                samplePPIs = json.loads(f.read())
            
            for i in range(len(samplePPIs)):
                fName = 'CRA_tenTrial_{}_{}'.format(dataset, i)
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(samplePPIs[i], "CRA", scoreArgsDict={}, coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def CRA_IntAct_spoke_tenTrial():
        for dataset in ['IntAct_spoke']:
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f:
                samplePPIs = json.loads(f.read())
            
            for i in range(len(samplePPIs)):
                fName = 'CRA_tenTrial_{}_{}'.format(dataset, i)
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(samplePPIs[i], "CRA", scoreArgsDict={}, coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

    @staticmethod
    def CAR_tenTrial():
        for dataset in ['bioGRID', 'STRING', 'MINT', 'IntAct_spoke']:
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f:
                samplePPIs = json.loads(f.read())
            
            for i in range(len(samplePPIs)):
                fName = 'CAR_tenTrial_{}_{}'.format(dataset, i)
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(samplePPIs[i], "CAR", scoreArgsDict={}, coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)
    
    @staticmethod
    def L3Raw_tenTrial():
        for dataset in ['bioGRID', 'STRING', 'MINT', 'IntAct_spoke']:
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f:
                samplePPIs = json.loads(f.read())
            
            for i in range(len(samplePPIs)):
                fName = 'L3Raw_tenTrial_{}_{}'.format(dataset, i)
                startTime = time.time()
                predPPI, predScore = ppiLPred.multiCore_PPILinkPred(samplePPIs[i], "L3Raw", scoreArgsDict={}, coreNo=14, logging=True)
                runTime = time.time()-startTime
                helper.write_runTime(fName, runTime)
                helper.write_resultData(predPPI, predScore, fName)

class dataTrimming:
    @staticmethod
    def load_tenTrial():
        tags = ['commonNeighbor']
        fNames = ['{}_tenTrial_bioGRID'.format(fName) for fName in tags]
        datasetClass = ['bioGRID' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

    @staticmethod
    def load_L3uvJoin_tenTrial():
        tags = ['L3uvJoin']
        fNames = ['{}_tenTrial_bioGRID'.format(fName) for fName in tags]
        datasetClass = ['bioGRID' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

    @staticmethod
    def load_ExactL3uvJoin_tenTrial():
        tags = ['xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        fNames = ['{}_tenTrial_bioGRID'.format(fName) for fName in tags]
        datasetClass = ['bioGRID' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)
    
    @staticmethod
    def load_MINT_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        fNames = ['{}_tenTrial_MINT'.format(fName) for fName in tags]
        datasetClass = ['MINT' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

    @staticmethod
    def load_IntAct_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        fNames = ['{}_tenTrial_IntAct'.format(fName) for fName in tags]
        datasetClass = ['IntAct' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)
    
    @staticmethod
    def load_IntAct_spoke_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        fNames = ['{}_tenTrial_IntAct_spoke'.format(fName) for fName in tags]
        datasetClass = ['IntAct_spoke' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

    @staticmethod
    def load_STRING_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        fNames = ['{}_tenTrial_STRING'.format(fName) for fName in tags]
        datasetClass = ['STRING' for tag in tags]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

    @staticmethod
    def load_CRA_tenTrial():
        datasetClass = ['bioGRID', 'STRING', 'MINT']
        fNames = ['CRA_tenTrial_{}'.format(ds) for ds in datasetClass]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)
    
    @staticmethod
    def load_CRA_IntAct_spoke_tenTrial():
        datasetClass = ['IntAct_spoke']
        fNames = ['CRA_tenTrial_{}'.format(ds) for ds in datasetClass]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

    @staticmethod
    def load_CAR_tenTrial():
        datasetClass = ['bioGRID', 'STRING', 'MINT', 'IntAct_spoke']
        fNames = ['CAR_tenTrial_{}'.format(ds) for ds in datasetClass]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)
    
    @staticmethod
    def load_L3Raw_tenTrial():
        datasetClass = ['bioGRID', 'STRING', 'MINT', 'IntAct_spoke']
        fNames = ['L3Raw_tenTrial_{}'.format(ds) for ds in datasetClass]
        helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

class PRCurveGen:
    @staticmethod
    def tenTrial_PRCurve():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        # split the processing to avoid lack of memory
        tags = ['commonNeighbor']

        # 40 for bioGRID
        print("bioGRID PR starts")
        fNames = ['{}_tenTrial_bioGRID'.format(fName) for fName in tags]
        bioGRID_samplePPIs = []
        with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [bioGRID_samplePPIs[i] for tag in tags for i in range(10)]
            , ['bioGRID' for tag in tags for i in range(10)], coreNo=13)
        print("bioGRID PR completed")

    @staticmethod
    def GGI_PRCurve():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['commonNeighbor']

        # 40 for bioGRID
        for tag in tags:
            fName = '{}_tenTrial_bioGRID'.format(tag)
            print('bioGRID '+fName+' started')
            bioGRID_samplePPIs = []
            with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore(["{}_GGI_{}".format(fName, i) for i in range(10)]
            , [predPPIs[fName][i] for i in range(10)], [bioGRID_samplePPIs[i] for i in range(10)]
            , ['bioGRID' for i in range(10)], coreNo=10, isGGI=True, logging=True)
            print('bioGRID '+fName+' completed')

    @staticmethod
    def L3uvJoin_tenTrial_PRCurve():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        # split the processing to avoid lack of memory
        tags = ['L3uvJoin']

        # 40 for bioGRID
        print("bioGRID PR starts")
        fNames = ['{}_tenTrial_bioGRID'.format(fName) for fName in tags]
        bioGRID_samplePPIs = []
        with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [bioGRID_samplePPIs[i] for tag in tags for i in range(10)]
            , ['bioGRID' for tag in tags for i in range(10)], coreNo=13)
        print("bioGRID PR completed")

    @staticmethod
    def ExactL3uvJoin_tenTrial_PRCurve():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']

        # 40 for bioGRID
        print("bioGRID PR starts")
        fNames = ['{}_tenTrial_bioGRID'.format(fName) for fName in tags]
        bioGRID_samplePPIs = []
        with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [bioGRID_samplePPIs[i] for tag in tags for i in range(10)]
            , ['bioGRID' for tag in tags for i in range(10)], coreNo=13)
        print("bioGRID PR completed")

    @staticmethod
    def uvJoin_GGIPR():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']

        # 40 for bioGRID
        for tag in tags:
            fName = '{}_tenTrial_bioGRID'.format(tag)
            print('bioGRID '+fName+' started')
            bioGRID_samplePPIs = []
            with open("./genData/bioGRID_sampledPPIs.json", "r") as f: bioGRID_samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore(["{}_GGI_{}".format(fName, i) for i in range(10)]
            , [predPPIs[fName][i] for i in range(10)], [bioGRID_samplePPIs[i] for i in range(10)]
            , ['bioGRID' for i in range(10)], coreNo=10, isGGI=True, logging=True)
            print('bioGRID '+fName+' completed')

    @staticmethod
    def MINT_PPIPR():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']

        print("MINT PR starts")
        fNames = ['{}_tenTrial_MINT'.format(fName) for fName in tags]
        MINT_samplePPIs = []
        with open("./genData/MINT_sampledPPIs.json", "r") as f: MINT_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [MINT_samplePPIs[i] for tag in tags for i in range(10)]
            , ['MINT' for tag in tags for i in range(10)], coreNo=13, logging=True)
        print("MINT PR completed")

    @staticmethod
    def IntAct_PPIPR():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']

        print("IntAct PR starts")
        fNames = ['{}_tenTrial_IntAct'.format(fName) for fName in tags]
        IntAct_samplePPIs = []
        with open("./genData/IntAct_sampledPPIs.json", "r") as f: IntAct_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [IntAct_samplePPIs[i] for tag in tags for i in range(10)]
            , ['IntAct' for tag in tags for i in range(10)], coreNo=13, logging=True)
        print("IntAct PR completed")
    
    @staticmethod
    def IntAct_spoke_PPIPR():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']

        print("IntAct spoke PR starts")
        fNames = ['{}_tenTrial_IntAct_spoke'.format(fName) for fName in tags]
        IntAct_samplePPIs = []
        with open("./genData/IntAct_spoke_sampledPPIs.json", "r") as f: IntAct_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [IntAct_samplePPIs[i] for tag in tags for i in range(10)]
            , ['IntAct_spoke' for tag in tags for i in range(10)], coreNo=13, logging=True)
        print("IntAct spoke PR completed")
    
    @staticmethod
    def STRING_PPIGGIPR():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']

        print("STRING PR starts")
        fNames = ['{}_tenTrial_STRING'.format(fName) for fName in tags]
        STRING_samplePPIs = []
        with open("./genData/STRING_sampledPPIs.json", "r") as f: STRING_samplePPIs = json.loads(f.read())
        helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
            , [predPPIs[fName][i] for fName in fNames for i in range(10)]
            , [STRING_samplePPIs[i] for tag in tags for i in range(10)]
            , ['STRING' for tag in tags for i in range(10)], coreNo=13, logging=True)
        print("STRING PR completed")

        for tag in tags:
            fName = '{}_tenTrial_STRING'.format(tag)
            print('STRING '+fName+' started')
            STRING_samplePPIs = []
            with open("./genData/STRING_sampledPPIs.json", "r") as f: STRING_samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore(["{}_GGI_{}".format(fName, i) for i in range(10)]
            , [predPPIs[fName][i] for i in range(10)], [STRING_samplePPIs[i] for i in range(10)]
            , ['STRING' for i in range(10)], coreNo=10, isGGI=True, logging=True)
            print('STRING '+fName+' completed')

    @staticmethod
    def CRA_PPIGGIPR():
        predPPIs  = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        print("CRA PR starts")
        for ds in ['bioGRID', 'STRING', 'MINT']:
            fName = 'CRA_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, logging=True)
        print("CRA PPI ends")

        print("CRA GGI starts")
        for ds in ['bioGRID', 'STRING']:
            fName = 'CRA_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_GGI_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, isGGI=True, logging=True)
        print("CRA GGI ends")
    
    @staticmethod
    def CRA_IntAct_spoke_PPIPR():
        predPPIs  = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        print("CRA PR starts")
        for ds in ['IntAct_spoke']:
            fName = 'CRA_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, logging=True)
        print("CRA PPI ends")

    @staticmethod
    def CAR_PPIGGIPR():
        predPPIs  = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        print("CAR PR starts")
        for ds in ['bioGRID', 'STRING', 'MINT', 'IntAct_spoke']:
            fName = 'CAR_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, logging=True)
        print("CAR PPI ends")

        print("CAR GGI starts")
        for ds in ['bioGRID', 'STRING']:
            fName = 'CAR_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_GGI_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, isGGI=True, logging=True)
        print("CAR GGI ends")
    
    @staticmethod
    def L3Raw_PPIGGIPR():
        predPPIs  = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        print("L3Raw PR starts")
        for ds in ['bioGRID', 'STRING', 'MINT', 'IntAct_spoke']:
            fName = 'L3Raw_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, logging=True)
        print("L3Raw PPI ends")

        print("L3Raw GGI starts")
        for ds in ['bioGRID', 'STRING']:
            fName = 'L3Raw_tenTrial_{}'.format(ds)
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(ds), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_GGI_{}".format(i) for i in range(10)]
                , predPPIs[fName], samplePPIs, [ds for i in range(10)], coreNo=10, isGGI=True, logging=True)
        print("L3Raw GGI ends")

class GoSemSIm_DataPrep:
    @staticmethod
    def tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['commonNeighbor']
        for dataset in ['bioGRID', 'STRING']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

    @staticmethod
    def uvJoin_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['bioGRID', 'STRING']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

    @staticmethod
    def IntAct_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['IntAct']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")
    
    @staticmethod
    def MINT_IntAct_spoke_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['MINT', 'IntAct_spoke']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")
    
    @staticmethod
    def STRING_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['STRING']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

    @staticmethod
    def CRA_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['CRA']
        for dataset in ['bioGRID', 'STRING', 'IntAct', 'MINT']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

    @staticmethod
    def CRA_IntAct_spoke_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['CRA']
        for dataset in ['IntAct_spoke']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

    @staticmethod
    def CAR_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['CAR']
        for dataset in ['bioGRID', 'STRING', 'IntAct_spoke', 'MINT']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")
    
    @staticmethod
    def L3Raw_tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['L3Raw']
        for dataset in ['bioGRID', 'STRING', 'IntAct_spoke', 'MINT']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = predPPIs[tag]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

class GOSemSim_Finalize:
    # after R script is ran on the prepared data in GOSemSimPrepData
    @staticmethod
    def tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI.json'):
            with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['commonNeighbor']
        for dataset in ['bioGRID', 'STRING']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

    @staticmethod
    def uvJoin_tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI.json'):
            with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['bioGRID', 'STRING']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

    @staticmethod
    def IntAct_tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI.json'):
            with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['IntAct']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))
    
    @staticmethod
    def MINT_IntAct_spoke_tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI.json'):
            with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['MINT', 'IntAct_spoke']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))
    
    @staticmethod
    def STRING_tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI.json'):
            with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin']
        for dataset in ['STRING']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))
    
    @staticmethod
    def CRA_tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI.json'):
            with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['CRA']
        for dataset in ['bioGRID', 'STRING', 'IntAct_spoke', 'MINT']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore.json', 'w') as f: f.write(json.dumps(GOScore))

if __name__ == "__main__":
    pass

    # initDataGen.CN_tenTrial()
    # initDataGen.L3_uvJoin_tenTrial()
    # initDataGen.ExactL3_uvJoin_tenTrial()
    # dataTrimming.load_tenTrial()
    # dataTrimming.load_L3uvJoin_tenTrial()
    # dataTrimming.load_ExactL3uvJoin_tenTrial()
    # PRCurveGen.tenTrial_PRCurve()
    # PRCurveGen.GGI_PRCurve()
    # PRCurveGen.L3uvJoin_tenTrial_PRCurve()
    # PRCurveGen.ExactL3uvJoin_tenTrial_PRCurve()
    # PRCurveGen.uvJoin_GGIPR()
    # GoSemSIm_DataPrep.tenTrial_to_csvData()
    # GoSemSIm_DataPrep.uvJoin_tenTrial_to_csvData()
    # GOSemSim_Finalize.tenTrial_res_to_json()
    # GOSemSim_Finalize.uvJoin_tenTrial_res_to_json()

    # initDataGen.IntAct_tenTrial()
    # dataTrimming.load_IntAct_tenTrial()
    # PRCurveGen.IntAct_PPIPR()
    # GoSemSIm_DataPrep.IntAct_tenTrial_to_csvData()
    # GOSemSim_Finalize.IntAct_tenTrial_res_to_json()

    # initDataGen.IntAct_spoke_tenTrial()
    # initDataGen.MINT_tenTrial()
    # dataTrimming.load_IntAct_spoke_tenTrial()
    # dataTrimming.load_MINT_tenTrial()
    # PRCurveGen.IntAct_spoke_PPIPR()
    # PRCurveGen.MINT_PPIPR()
    # GoSemSIm_DataPrep.MINT_IntAct_spoke_tenTrial_to_csvData()

    # initDataGen.STRING_tenTrial()
    # dataTrimming.load_STRING_tenTrial()
    # PRCurveGen.STRING_PPIGGIPR()
    # GoSemSIm_DataPrep.STRING_tenTrial_to_csvData() # annotate each R script ran

    # GOSemSim_Finalize.MINT_IntAct_spoke_tenTrial_res_to_json()
    # GOSemSim_Finalize.STRING_tenTrial_res_to_json()

    # initDataGen.CRA_tenTrial()
    # dataTrimming.load_CRA_tenTrial()
    # PRCurveGen.CRA_PPIGGIPR()
    # GoSemSIm_DataPrep.CRA_tenTrial_to_csvData()
    # initDataGen.CRA_IntAct_spoke_tenTrial()
    # dataTrimming.load_CRA_IntAct_spoke_tenTrial()
    # PRCurveGen.CRA_IntAct_spoke_PPIPR()
    # GoSemSIm_DataPrep.CRA_IntAct_spoke_tenTrial_to_csvData()

    # initDataGen.CAR_tenTrial()
    # dataTrimming.load_CAR_tenTrial()
    # PRCurveGen.CAR_PPIGGIPR()
    # GoSemSIm_DataPrep.CAR_tenTrial_to_csvData()

    # initDataGen.L3Raw_tenTrial()
    # dataTrimming.load_L3Raw_tenTrial()
    # PRCurveGen.L3Raw_PPIGGIPR()
    # GoSemSIm_DataPrep.L3Raw_tenTrial_to_csvData()

    # GOSemSim_Finalize.CRA_tenTrial_res_to_json()