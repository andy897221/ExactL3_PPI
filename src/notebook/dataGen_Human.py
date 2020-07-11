import sys
sys.path.append('..')
import json
import os
from itertools import combinations
import onPPILinkPred_HPC as ppiLPred
import bioGRID as bg
import HuRI
import IntAct
import STRING as string
import numpy as np
import traversalHelper as tr
import pandas as pd
import domain as dm
from collections import defaultdict
import time
import MINT
import helper as hr

class helper:
    @staticmethod
    def write_runTime(tag, runTime):
        fName = './resultData/homo/runTime.json'
        if not os.path.exists(fName):
            with open(fName, 'w') as f: f.write(json.dumps({}))
        
        runTimeD = {}
        with open(fName, 'r') as f: runTimeD = json.loads(f.read())
        runTimeD[tag] = runTime
        with open(fName, 'w') as f: f.write(json.dumps(runTimeD))

    @staticmethod
    def trim_multiple_ppi_result(fNames, datasetClass, trialSize):
        bioGRID_homo = int(len([*bg.parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.187.tab2.txt'
            , wFile_GGI='./data/parsed/BioGRID_homo_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_homo_PPI.pkl', root="../")][1].index)*0.5)

        STRING_homo = int(len([*string.parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
            , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
            , uniProtMap='./data/UniProt/uniprot-taxonomy_9606_STRING.tab', root='../'
            , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')][1].index)*0.5)
        
        MINT_homo = int(len([*MINT.parse_MINT(ppiFile='./data/MINT/species human', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
                wFile_GGI='./data/parsed/MINT_homo_GGI.pkl', wFile_PPI='./data/parsed/MINT_homo_PPI.pkl', root="../")][1].index)*0.5)

        trimNum = {'HuRI': int(len(HuRI.parse_HuRI(root='../').index)*0.5)
        , "bioGRID_homo": bioGRID_homo, "STRING_homo": STRING_homo, "MINT_homo": MINT_homo}

        if not os.path.exists('./resultData/homo/trimmed_predPPIs.json'):
            with open('./resultData/homo/trimmed_predPPIs.json', 'w') as f: pass
            with open('./resultData/homo/trimmed_predScores.json', 'w') as f: pass
        
        for i in range(len(fNames)):
            predPPI, predScore = [], []
            for j in range(trialSize):
                with open("./resultData/homo/{}_{}_PPI.json".format(fNames[i], j), 'r') as f:
                    predPPI.append(json.loads(f.read())[0:trimNum[datasetClass[i]]])
                with open("./resultData/homo/{}_{}_score.json".format(fNames[i], j), 'r') as f:
                    predScore.append(json.loads(f.read())[0:trimNum[datasetClass[i]]])
            predPPIs, predScores = {fNames[i]: predPPI}, {fNames[i]: predScore}
            with open('./resultData/homo/trimmed_predPPIs.json', 'a+') as f: f.write(json.dumps(predPPIs)+"\n")
            with open('./resultData/homo/trimmed_predScores.json', 'a+') as f: f.write(json.dumps(predScores)+"\n")

    @staticmethod
    def write_resultData(predPPI, predScore, fName):
        with open('./resultData/homo/{}_PPI.json'.format(fName), 'w') as f:
            f.write(json.dumps(predPPI))
        with open('./resultData/homo/{}_score.json'.format(fName), 'w') as f:
            f.write(json.dumps(predScore))

    @staticmethod
    def append_precRecMap_multiCore(fNames, predPPI, samplePPI, datasetClass, coreNo, isGGI=False, logging=False):
        if isGGI: i = 0
        else: i = 1

        bioGRID_homo = [list(ppi) for ppi in np.asarray([*bg.parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.187.tab2.txt'
            , wFile_GGI='./data/parsed/BioGRID_homo_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_homo_PPI.pkl', root="../")][i][['nodeA', 'nodeB']])]

        STRING_homo = [list(ppi) for ppi in np.asarray([*string.parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
            , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
            , uniProtMap='./data/UniProt/uniprot-taxonomy_9606_STRING.tab', root='../'
            , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')][i][['nodeA', 'nodeB']])]
        
        MINT_homo = [list(ppi) for ppi in np.asarray([*MINT.parse_MINT(ppiFile='./data/MINT/species human', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
                wFile_GGI='./data/parsed/MINT_homo_GGI.pkl', wFile_PPI='./data/parsed/MINT_homo_PPI.pkl', root="../")][i][['nodeA', 'nodeB']])]


        fullPPISet = {'HuRI': [list(ppi) for ppi in np.asarray(HuRI.parse_HuRI(root='../')[['nodeA', 'nodeB']])]
        , 'bioGRID_homo': bioGRID_homo, 'STRING_homo': STRING_homo, 'MINT_homo': MINT_homo}

        fullPrecRecMap = {}
        if not os.path.exists("./resultData/PRCurveMap_homo.json"):
            with open("./resultData/PRCurveMap_homo.json", "w") as f: f.write(json.dumps(fullPrecRecMap))

        precRecMap = ppiLPred.precRecMap_multiCore(fNames, predPPI, samplePPI, [fullPPISet[i] for i in datasetClass], coreNo, logging)
        with open('./resultData/PRCurveMap_homo.json', 'r') as f: fullPrecRecMap = json.loads(f.read())
        fullPrecRecMap.update(precRecMap)
        with open('./resultData/PRCurveMap_homo.json', 'w') as f: f.write(json.dumps(fullPrecRecMap))

    @staticmethod
    def uniprot_map(file="../data/UniProt/uniprot-taxonomy_9606.tab"):
        df = pd.read_csv(file, sep='\t')
        df.dropna(inplace=True)
        genes = df["Gene names  (primary )"]
        entry = df["Entry"]
        geneToEntry = dict(zip(list(genes), list(entry)))
        return geneToEntry

class dataGen:
    @staticmethod
    def homo_HPC():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
        methods = ['commonNeighbor', 'L3uvJoin', 'interStr', 'CRA', 'Sim', 'CH2_L3']
        scoreArgDicts = [{}, {}, {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}, {}, {}, {}]

        datasets = ["STRING_homo", "bioGRID_homo"]

        for dataset in datasets:
            for t in range(len(tags)):
                for i in range(10):
                    samplePPI = []
                    with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f: samplePPI = json.loads(f.read())
                    samplePPI = samplePPI[i]

                    fName = '{}_tenTrial_{}_{}'.format(tags[t], dataset, i)
                    startTime = time.time()
                    predPPI, predScore, rank = ppiLPred.multiCore_PPILinkPred_HPC(samplePPI, methods[t], scoreArgsDict=scoreArgDicts[t], coreNo=14, logging=False)
                    runTime = time.time()-startTime
                    if rank == 0: helper.write_runTime(fName, runTime)
                    helper.write_resultData(predPPI, predScore, fName+"_c{}".format(rank))

    @staticmethod
    def MINT_HuRI_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']

        for dataset in ['HuRI', 'MINT_homo']:
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
                    
    @staticmethod
    def candidateSz_all():
        datasets = ["bioGRID_homo", "STRING_homo"]
        for dataset in datasets:
            for i in range(10):
                samplePPI = []
                with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f: samplePPI = json.loads(f.read())
                samplePPI = samplePPI[i]

                count, rank = ppiLPred.multiCore_PPILinkPred_count(samplePPI, coreNo=24)
                fName = 'tenTrial_{}_{}_c{}'.format(dataset, i, rank)
                with open('./resultData/homo/{}_candidateNum.txt'.format(fName), 'w') as f:
                    f.write(str(count))

                    
class dataTrimming:
    @staticmethod
    def STRING_homo_HPC_combine():
        dataset_len = int(len([*string.parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
            , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
            , uniProtMap='./data/UniProt/uniprot-taxonomy_9606_STRING.tab', root='../'
            , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')][1].index)*0.5)

        dataset = "STRING"
        datasetName = "STRING_homo"
        for folder in os.listdir("./resultData/homo/"):
            if dataset in folder:
                for trial in range(0, 10):
                    print("folder: {}, trial: {}".format(folder, trial))
                    filenames = os.listdir("./resultData/homo/"+folder)
                    ppiFiles = [i for i in filenames if "PPI" in i]
                    curTrialFiles = [i for i in ppiFiles if i.split("_")[-3] == str(trial)]

                    topPPIs, topScores = [], []
                    for i in curTrialFiles:
                        with open("./resultData/homo/{}/{}".format(folder, i), "r") as f:
                            topPPIs += json.loads(f.read())
                        with open("./resultData/homo/{}/{}".format(folder, "_".join(i.split("_")[:-1])+"_score.json"), "r") as f:
                            topScores += json.loads(f.read())
                        # sort
                        topPPIs, topScores = hr.sort_key_val(topPPIs, topScores)
                        topPPIs = topPPIs[0:dataset_len]
                        topScores = topScores[0:dataset_len]

                    with open("./resultData/homo/{}/{}_{}_trimmedPPIs.json".format(folder,curTrialFiles[0].split("_")[0], datasetName), "a+") as f:
                        f.write(json.dumps(topPPIs)+"\n")
                    with open("./resultData/homo/{}/{}_{}_trimmedScores.json".format(folder,curTrialFiles[0].split("_")[0], datasetName), "a+") as f:
                        f.write(json.dumps(topScores)+"\n")

    @staticmethod
    def group_STRING_homo():
        tagMap = {
            "CH2": "CH2_L3",
            "CN": "commonNeighbor",
            'ExactL3': "xyContrib_dualCN_uvJoin",
            "Sim": "Sim",
            "CRA": "CRA",
            "L3uvJoin": "L3uvJoin"
        }

        dataset = "STRING"
        datasetName = "STRING_homo"
        ppiFiles, scoreFiles = {}, {}
        for folder in os.listdir("./resultData/homo/"):
            if dataset in folder:
                filenames = os.listdir("./resultData/homo/"+folder)
                
                ppis = []
                with open(["./resultData/homo/{}/{}".format(folder, i) for i in filenames if "trimmedPPIs" in i][0], "r") as f:
                    for line in f.readlines(): ppis.append(json.loads(line))
                ppiFiles["{}_tenTrial_{}".format(tagMap[folder.split("_")[0]], datasetName)] = ppis
                
                scores = []
                with open(["./resultData/homo/{}/{}".format(folder, i) for i in filenames if "trimmedScores" in i][0], "r") as f:
                    for line in f.readlines(): scores.append(json.loads(line))
                scoreFiles["{}_tenTrial_{}".format(tagMap[folder.split("_")[0]], datasetName)] = scores

        with open("./resultData/homo/trimmed_predPPIs_homo.json", "a+") as f:
            f.write(json.dumps(ppiFiles)+"\n")
        with open("./resultData/homo/trimmed_predScores_homo.json", "a+") as f:
            f.write(json.dumps(scoreFiles)+"\n")

    @staticmethod
    def bioGRID_homo_HPC_combine():
        dataset_len = int(len([*bg.parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.187.tab2.txt'
            , wFile_GGI='./data/parsed/BioGRID_homo_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_homo_PPI.pkl', root="../")][1].index)*0.5)

        dataset = "bioGRID"
        datasetName = "bioGRID_homo"
        for folder in os.listdir("./resultData/homo/"):
            if dataset in folder:
                for trial in range(6, 10):
                    print("folder: {}, trial: {}".format(folder, trial))
                    filenames = os.listdir("./resultData/homo/"+folder)
                    ppiFiles = [i for i in filenames if "PPI" in i]
                    curTrialFiles = [i for i in ppiFiles if i.split("_")[-3] == str(trial)]

                    topPPIs, topScores = [], []
                    for i in curTrialFiles:
                        with open("./resultData/homo/{}/{}".format(folder, i), "r") as f:
                            topPPIs += json.loads(f.read())
                        with open("./resultData/homo/{}/{}".format(folder, "_".join(i.split("_")[:-1])+"_score.json"), "r") as f:
                            topScores += json.loads(f.read())
                        # sort
                        topPPIs, topScores = hr.sort_key_val(topPPIs, topScores)
                        topPPIs = topPPIs[0:dataset_len]
                        topScores = topScores[0:dataset_len]

                    with open("./resultData/homo/{}/{}_{}_trimmedPPIs.json".format(folder,curTrialFiles[0].split("_")[0], datasetName), "a+") as f:
                        f.write(json.dumps(topPPIs)+"\n")
                    with open("./resultData/homo/{}/{}_{}_trimmedScores.json".format(folder,curTrialFiles[0].split("_")[0], datasetName), "a+") as f:
                        f.write(json.dumps(topScores)+"\n")

    @staticmethod
    def group_bioGRID_homo():
        tagMap = {
            "CH2": "CH2_L3",
            "CN": "commonNeighbor",
            'ExactL3': "xyContrib_dualCN_uvJoin",
            "Sim": "Sim",
            "CRA": "CRA",
            "L3uvJoin": "L3uvJoin"
        }
            
        datasetNames = ['bioGRID_homo', 'STRING_homo']
        datasets = ['bioGRID', 'STRING']
        for k in range(datasets):
            dataset = datasets[k]
            datasetName = datasetNames[k]
            
            ppiFiles, scoreFiles = {}, {}
            for folder in os.listdir("./resultData/homo/"):
                if dataset in folder:
                    filenames = os.listdir("./resultData/homo/"+folder)
                    
                    ppis = []
                    with open(["./resultData/homo/{}/{}".format(folder, i) for i in filenames if "trimmedPPIs" in i][0], "r") as f:
                        for line in f.readlines(): ppis.append(json.loads(line))
                    ppiFiles["{}_tenTrial_{}".format(tagMap[folder.split("_")[0]], datasetName)] = ppis
                    
                    scores = []
                    with open(["./resultData/homo/{}/{}".format(folder, i) for i in filenames if "trimmedScores" in i][0], "r") as f:
                        for line in f.readlines(): scores.append(json.loads(line))
                    scoreFiles["{}_tenTrial_{}".format(tagMap[folder.split("_")[0]], datasetName)] = scores

            with open("./resultData/homo/trimmed_predPPIs_homo.json", "a+") as f:
                f.write(json.dumps(ppiFiles)+"\n")
            with open("./resultData/homo/trimmed_predScores_homo.json", "a+") as f:
                f.write(json.dumps(scoreFiles)+"\n")

    @staticmethod
    def load_tenTrial():
        tags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
        for dataset in ['HuRI', 'MINT_homo']:
            fNames = ['{}_tenTrial_{}'.format(fName, dataset) for fName in tags]
            datasetClass = [dataset for tag in tags]
            helper.trim_multiple_ppi_result(fNames, datasetClass, 10)
                
class PRCurveGen:
    @staticmethod
    def homo_PPIPR():
        predPPIs = {}
        with open('./resultData/homo/trimmed_predPPIs_homo.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))

        tags_each = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3']
        datasets = ['bioGRID_homo', 'STRING_homo', 'MINT_homo', 'HuRI']
        
        for datset in datasets:
            for tags in tags_each:
                print("{} PR starts".format(dataset))
                tags = [tags]
                fNames = ['{}_tenTrial_{}'.format(fName, dataset) for fName in tags]
                samplePPIs = []
                with open("./genData/{}_sampledPPIs.json".format(dataset), "r") as f: samplePPIs = json.loads(f.read())
                helper.append_precRecMap_multiCore([fName+"_{}".format(i) for fName in fNames for i in range(10)]
                    , [predPPIs[fName][i] for fName in fNames for i in range(10)]
                    , [samplePPIs[i] for tag in tags for i in range(10)]
                    , [dataset for tag in tags for i in range(10)], coreNo=8, logging=True)
            print("{} PR completed".format(dataset))

class GOSemSim_DataPrep:
    @staticmethod
    def HuRI_tenTrial_to_csvData():
        predPPIs = {}
        HuRI_df = HuRI.parse_HuRI(root="../")
        dataset_len = int(int(len(HuRI_df.index)*0.5)*0.1)
        with open('./resultData/homo/trimmed_predPPIs_homo.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        geneToEntry = helper.uniprot_map()
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'CH2_L3', 'Sim']
        for dataset in ['HuRI']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = [[[geneToEntry[g[0]], geneToEntry[g[1]]] for g in j[:dataset_len] if g[0] in geneToEntry and g[1] in geneToEntry] for j in predPPIs[tag]]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")
    
    @staticmethod
    def MINT_homo_tenTrial_to_csvData():
        predPPIs = {}
        ggi_df, ppi_df = MINT.parse_MINT(ppiFile='./data/MINT/species human', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
            wFile_GGI='./data/parsed/MINT_homo_GGI.pkl', wFile_PPI='./data/parsed/MINT_homo_PPI.pkl', root='../')
        dataset_len = int(int(len(ppi_df.index)*0.5)*0.1)
        with open('./resultData/homo/trimmed_predPPIs_homo.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        geneToEntry = helper.uniprot_map()
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'CH2_L3', 'Sim']
        for dataset in ["MINT_homo"]:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = [[[geneToEntry[g[0]], geneToEntry[g[1]]] for g in j[:dataset_len] if g[0] in geneToEntry and g[1] in geneToEntry] for j in predPPIs[tag]]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

    @staticmethod
    def STRING_homo_tenTrial_to_csvData():
        predPPIs = {}
        ggi_df, ppi_df = string.parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
        , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
        , uniProtMap='./data/UniProt/uniprot-taxonomy_9606_STRING.tab', root='../'
        , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')
        dataset_len = int(int(len(ppi_df.index)*0.5)*0.1)
        with open('./resultData/homo/trimmed_predPPIs_homo.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
        geneToEntry = helper.uniprot_map()
        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'CH2_L3', 'Sim']
        for dataset in ["STRING_homo"]:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = [[[geneToEntry[g[0]], geneToEntry[g[1]]] for g in j[:dataset_len] if g[0] in geneToEntry and g[1] in geneToEntry] for j in predPPIs[tag]]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")
    
    @staticmethod
    def bioGRID_homo_tenTrial_to_csvData():
        bioGRID_GGI_df, bioGRID_PPI_df = bg.parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.187.tab2.txt'
        , wFile_GGI='./data/parsed/BioGRID_homo_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_homo_PPI.pkl', root="../")
        dataset_len = int(int(len(bioGRID_PPI_df.index)*0.5)*0.1)
        geneToEntry = helper.uniprot_map()
        predPPIs = {}
        with open('./resultData/homo/trimmed_predPPIs_homo.json', 'r') as f:
            for line in f.readlines(): predPPIs.update(json.loads(line))
            
        baseTags = ['commonNeighbor', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3', 'L3uvJoin']
        for dataset in ["bioGRID_homo"]:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for i in range(len(tags)):
                tag = tags[i]
                PPIsList = [[[geneToEntry[g[0]], geneToEntry[g[1]]] for g in j[:dataset_len] if g[0] in geneToEntry and g[1] in geneToEntry] for j in predPPIs[tag]]
                for j in range(len(PPIsList)):
                    with open('./GoSemSimPrepData/{}_{}.csv'.format(tag, j), 'w') as f:
                        f.write("\n".join(['nodeA\tnodeB']+["\t".join(ppi) for ppi in PPIsList[j]])+"\n")

class GOSemSim_Finalize:
    @staticmethod
    def homo_tenTrial_res_to_json():
        PPI, GOScore = {}, {}
        if not os.path.exists('./GOSemSimFinalizedData/GOPPI_homo.json'):
            with open('./GOSemSimFinalizedData/GOPPI_homo.json', 'w') as f: f.write(json.dumps(PPI))
            with open('./GOSemSimFinalizedData/GOScore_homo.json', 'w') as f: f.write(json.dumps(GOScore))

        with open('./GOSemSimFinalizedData/GOPPI_homo.json', 'r') as f: PPI = json.loads(f.read())
        with open('./GOSemSimFinalizedData/GOScore_homo.json', 'r') as f: GOScore = json.loads(f.read())

        baseTags = ['commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'CH2_L3', 'Sim']
        for dataset in ['STRING_homo', 'HuRI', 'bioGRID_homo', 'MINT_homo']:
            tags = ["{}_tenTrial_{}".format(baseTag, dataset) for baseTag in baseTags]
            for tag in tags:
                PPI[tag], GOScore[tag] = [], []
                for i in range(0, 10):
                    fName = './GOSemSimResData/{}_{}.csv'.format(tag, i)
                    df = pd.read_csv(fName, sep='\t')
                    PPI[tag].append([list(ppi) for ppi in np.asarray(df[['nodeA', 'nodeB']])])
                    GOScore[tag].append([float(score) for score in df['score']])
        
        with open('./GOSemSimFinalizedData/GOPPI_homo.json', 'w') as f: f.write(json.dumps(PPI))
        with open('./GOSemSimFinalizedData/GOScore_homo.json', 'w') as f: f.write(json.dumps(GOScore))

if __name__ == "__main__":
    pass
    # run HuRI and MINT_homo
    # dataGen.MINT_HuRI_tenTrial()
    # dataTrimming.load_tenTrial()
    # GOSemSim_DataPrep.HuRI_tenTrial_to_csvData()
    # GOSemSim_DataPrep.MINT_homo_tenTrial_to_csvData()
    
    # run STRING and bioGRID on high cluster computer using mpi4py
    # dataGen.homo_HPC()
    # dataTrimming.STRING_homo_HPC_combine()
    # dataTrimming.bioGRID_homo_HPC_combine()
    # dataTrimming.group_bioGRID_homo()
    # dataTrimming.group_STRING_homo()
    # PRCurveGen.homo_PPIPR()
    # GOSemSim_DataPrep.bioGRID_homo_tenTrial_to_csvData()
    # GOSemSim_DataPrep.STRING_homo_tenTrial_to_csvData()
    
    # after ran GOSemSim R Scripts
    # GOSemSim_Finalize.homo_tenTrial_res_to_json()
    
    # gen PRCurve
    # PRCurveGen.homo_PPIPR()
    
    # evaluate candidate edges number
    # dataGen.candidateSz_all()