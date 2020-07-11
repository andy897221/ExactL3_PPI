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
    def tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']

        for dataset in ['bioGRID', 'STRING', 'IntAct_spoke', 'MINT']:
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f: samplePPIs = json.loads(f.read())
            for t in range(len(tags)):
                for i in range(len(samplePPIs)):
                    samplePPI = samplePPIs[i]
                    startTime = time.time()
                    predPPI, predScore = ppiLPred.multiCore_PPILinkPred(samplePPI, tags[t], scoreArgsDict={}, coreNo=14, logging=True)
                    runTime = time.time()-startTime
                    helper.write_runTime('{}_tenTrial_{}_{}'.format(tags[t], dataset, i), runTime)
                    helper.write_resultData(predPPI, predScore, '{}_tenTrial_{}_{}'.format(tags[t], dataset, i))

class dataTrimming:
    @staticmethod
    def load_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
        for dataset in ['bioGRID', 'STRING', 'IntAct_spoke', 'MINT']:
            fNames = ['{}_tenTrial_{}'.format(fName, dataset) for fName in tags]
            datasetClass = [dataset for tag in tags]
            helper.trim_multiple_ppi_result(fNames, datasetClass, 10)

class PRCurveGen:
    @staticmethod
    def tenTrial_PRCurve():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
        for dataset in ['bioGRID', 'STRING', 'IntAct_spoke', 'MINT']:
            print("{} PR starts".format(dataset))
            fNames = ['{}_tenTrial_{}'.format(fName, dataset) for fName in tags]
            samplePPIs = []
            with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f: samplePPIs = json.loads(f.read())
            helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
                , [predPPIs[fName][i] for fName in fNames for i in range(10)]
                , [samplePPIs[i] for tag in tags for i in range(10)]
                , [dataset for tag in tags for i in range(10)], coreNo=13)
            print("{} PR completed".format(dataset))

class GoSemSIm_DataPrep:
    @staticmethod
    def tenTrial_to_csvData():
        predPPIs = {}
        with open('./resultData/trimmed_predPPIs.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
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

        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
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
    # workflow:
    # run initDataGen.tenTrial() to gen tenTrial of n dataset with m algorithms
    # run dataTrimming.tenTrial() to parse those data generated into another file with reduced size (sampled)
    # run PRCurveGen.tenTrial_PRCruve() to gen Precision-Recall Curve based on the ranked PPI data against the original database
    # run GOSemSIm_DataPrep.tenTrial_to_csvData() to gen csv file from the ranked PPI into a format that can be processed by the R script to generate GOSemSim scores of those PPIs
    # run GOSemSim_Finalize.tenTrial_res_to_json() to parse those generated GOSemSim scores from the R script (ran seperately) into json format that can be processed by my python notebook
    
    # initDataGen.tenTrial()
    # dataTrimming.load_tenTrial()
    # PRCurveGen.tenTrial_PRCurve()
    # GOSemSIm_DataPrep.tenTrial_to_csvData()
    # GOSemSim_Finalize.tenTrial_res_to_json()