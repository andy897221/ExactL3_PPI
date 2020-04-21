import os
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import pickle
from itertools import combinations
from collections import defaultdict
import helper

# entryset -> entry -> 
# interactionList -> participantList -> participant  -> names -> shortLabel
#                                                    -> interactorRef (text) (id of interactor)
#                                                    -> experimentalRoleList -> experimentalRole -> names -> shortLabel (text)
#                 -> interactionType -> names -> shortLabel -> (text e.g. physical association)
#                                             -> fullname -> (text)
# interactorList -> interactor (<id=''>) -> names -> shortLabel -> (text (remove '_yeast'))

def parse_all_elem(db):
    elemList = []
    for elem in db.iter():
        elemList.append(elem.tag)
    elemList = list(set(elemList))
    return elemList

def get_prefix(db):
    return parse_all_elem(db)[0].split("}")[0]+"}"

def get_tag(prefix, tag):
    return prefix+tag

def parse_IntAct(folderName="./data/IntAct/yeast/", root="./"
    , wFile_GGI="./data/parsed/IntAct_GGI", wFile_PPI="./data/parsed/IntAct_PPI", spokeModel=False):

    if spokeModel: wFile_GGI, wFile_PPI = wFile_GGI+'_spoke.pkl', wFile_PPI+'_spoke.pkl'
    else: wFile_GGI, wFile_PPI = wFile_GGI+'.pkl', wFile_PPI+'.pkl'
    wFile_GGI, wFile_PPI, folderName = root+wFile_GGI, root+wFile_PPI, root+folderName

    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)
    
    PPIs, PPITypes = [], []
    ID_PPIName = {}
    for ppiF in os.listdir(folderName):
        db = ET.parse(folderName+"/"+ppiF).getroot()
        prefix = get_prefix(db)
        db = db.findall(get_tag(prefix, "entry"))[0]
        
        interactor = db.findall(get_tag(prefix, "interactorList"))[0]
        for node in interactor.findall(get_tag(prefix, "interactor")):
            nodeID = node.get("id")
            nodeName = node.findall(get_tag(prefix, "names"))[0].findall(get_tag(prefix, "shortLabel"))[0].text.upper()
            ID_PPIName[nodeID] = nodeName
        
        interaction = db.findall(get_tag(prefix, "interactionList"))[0]
        for ppis in interaction.findall(get_tag(prefix, "interaction")):
            ppiNodes = ppis.findall(get_tag(prefix, "participantList"))[0]
            ppiName = []
            role = []
            for protein in ppiNodes.findall(get_tag(prefix, "participant")):
                # print(ppiF, protein.get('id'))
                ppiName.append(protein.findall(get_tag(prefix, "interactorRef"))[0].text)
                role.append(protein.findall(get_tag(prefix, "experimentalRoleList"))[0].findall(get_tag(prefix, "experimentalRole"))[0].findall(get_tag(prefix, "names"))[0].findall(get_tag(prefix, "shortLabel"))[0].text.lower())
            ppiRoles = dict(zip(ppiName, role))
            if len(ppiName) < 2: continue # remove non binary PPI (mono PPI (self interact?))
            if ppiName[0] == ppiName[1]: continue # remove self-interact
            ppiName = sorted(ppiName)
            if ppiName in PPIs: continue # remove directional

            PPIType = ppis.findall(get_tag(prefix, "interactionType"))[0].findall(get_tag(prefix, "names"))[0].findall(get_tag(prefix, "shortLabel"))[0].text
            if spokeModel:
                modelMap = defaultdict(list)
                for p in ppiRoles: modelMap[ppiRoles[p]].append(p)
                ppiName = [[x,y] for x in list(modelMap['prey']) for y in list(modelMap['bait'])]

                # matrix model for any unspecified proteins
                remainingP = []
                for r in modelMap:
                    if r not in ['prey', 'bait']: remainingP += modelMap[r]
                ppiName += list(combinations(remainingP, 2))

                PPIType = [PPIType]*len(ppiName)
            else:
                ppiName = list(combinations(ppiName, 2))
                ppiName = [i for i in ppiName if i[0] != i[1]] # remove self-interact
                PPIType = [PPIType]*len(ppiName)
            PPIs, PPITypes = PPIs+ppiName, PPITypes+PPIType

    for j in range(len(PPIs)):
        PPIs[j] = np.asarray([ID_PPIName[i] for i in PPIs[j]])
    PPIs = np.transpose(np.asarray(PPIs))
    df = pd.DataFrame({'nodeA': PPIs[0], 'nodeB': PPIs[1], 'type':PPITypes})
    df = df[df.nodeA.str.contains("_YEAST", case=False)]
    df = df[df.nodeB.str.contains("_YEAST", case=False)]

    # remove _YEAST postfix
    df['nodeA'] = pd.Series(["".join(node.split("_YEAST")) for node in df['nodeA']]).values
    df['nodeB'] = pd.Series(["".join(node.split("_YEAST")) for node in df['nodeB']]).values
    df.reset_index(inplace=True, drop=True)

    # print(set(list(df['type'])))
    # {'physical association', 'methylation', 'association', 'colocalization', 'disulfide bond', 'phosphorylation', 'dephosphorylation', 'deubiquitination', 'enzymatic reaction', 'direct interaction'}

    ppi_df = df[df['type'].isin(["physical association", "association", "colocalization", "direct interaction"])].copy()
    ggi_df = df[~df['type'].isin(["physical association", "association", "colocalization", "direct interaction"])].copy()
    ppi_df.reset_index(inplace=True, drop=True)
    ggi_df.reset_index(inplace=True, drop=True)

    ppi_df.to_pickle(wFile_PPI)
    ggi_df.to_pickle(wFile_GGI)
    return ggi_df, ppi_df


if __name__ == "__main__":
    ggi_df, ppi_df = parse_IntAct()
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ggi_df, ppi_df = parse_IntAct(spokeModel=True)
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))