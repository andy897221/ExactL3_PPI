import pandas as pd
import bioGRID as bg
import traversalHelper as tr
import numpy as np
import os
from collections import defaultdict
from statistics import mean

def parse_uniProt_map(uniProtMapF):
    df = pd.read_csv(uniProtMapF, sep='\t')
    df.dropna(inplace=True)
    uniProtMapping = dict(zip([i.split(";")[0] for i in df['Cross-reference (STRING)']], list(df['Gene names  (primary )'])))
    return uniProtMapping

def parse_STRING(ppiFile='./data/STRING/4932.protein.links.v11.0.txt'
    , typeFile='./data/STRING/4932.protein.actions.v11.0.txt'
    , uniProtMap='./data/UniProt/uniprot-taxonomy_559292_STRING.tab', root='./'
    , wFile_GGI='./data/parsed/STRING_GGI.pkl', wFile_PPI='./data/parsed/STRING_PPI.pkl'):
    ppiFile, typeFile, wFile_GGI, wFile_PPI, uniProtMap = root+ppiFile, root+typeFile, root+wFile_GGI, root+wFile_PPI, root+uniProtMap

    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)

    # Sys name (used by STRING) => gene name (used by this project)
    reverseGeneMap = parse_uniProt_map(uniProtMap)

    df_STRING = pd.read_csv(ppiFile, sep=' ')
    dfType_STRING = pd.read_csv(typeFile, sep='\t')

    ppi, ppiScore = [], [] # build STRING ppi with gene name standard
    for i in range(len(df_STRING['protein1'])):
        p1, p2 = df_STRING['protein1'][i], df_STRING['protein2'][i]
        if p1 not in reverseGeneMap or p2 not in reverseGeneMap: continue
        ppi.append([reverseGeneMap[p1], reverseGeneMap[p2]])
        ppiScore.append(df_STRING['combined_score'][i])

    # make scores mean of conf scores of ppi of both direction
    sortedPPI = tr.Helper.list_to_pathStrs([sorted(i) for i in ppi])
    dualScoreMap = defaultdict(list)
    for i in range(len(sortedPPI)): dualScoreMap[sortedPPI[i]].append(ppiScore[i])

    ppiTypeMap = {} # parse STRING ppi type
    for i in range(len(dfType_STRING['item_id_a'])):
        p1, p2 = dfType_STRING['item_id_a'][i], dfType_STRING['item_id_b'][i]
        if p1 not in reverseGeneMap or p2 not in reverseGeneMap: continue
        if p1 == p2: continue
        ppiTypeMap[tr.Helper.br_to_pathStr([reverseGeneMap[p1], reverseGeneMap[p2]])] = dfType_STRING['mode'][i]
        ppiTypeMap[tr.Helper.br_to_pathStr([reverseGeneMap[p2], reverseGeneMap[p1]])] = dfType_STRING['mode'][i]

    # standard PPI DataFrame building
    ppiKeys = [ppiKey for ppiKey in list(dualScoreMap.keys()) if ppiKey in ppiTypeMap]
    ppi = list(np.transpose(np.asarray(tr.Helper.pathStrs_to_list(ppiKeys))))
    scores = [float(mean(dualScoreMap[ppiKey])/1000) for ppiKey in ppiKeys]
    ppiType = [ppiTypeMap[ppiKey] for ppiKey in ppiKeys]
    parsed_df = pd.DataFrame({'nodeA': ppi[0], 'nodeB': ppi[1], 'type': ppiType, 'score': scores})

    ggi_df = parsed_df[parsed_df['type'] != 'binding']
    ppi_df = parsed_df[parsed_df['type'] == 'binding']

    ggi_df.to_pickle(wFile_GGI)
    ppi_df.to_pickle(wFile_PPI)
    return ggi_df, ppi_df

if __name__ == "__main__":
    ggi_df, ppi_df = parse_STRING()
    # print(set(list(STRING_GGI['type'])))
    # {'catalysis', 'ptmod', 'inhibition', 'activation', 'reaction', 'expression'}
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))

    ggi_df, ppi_df = parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
    , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
    , uniProtMap='./data/UniProt/uniprot-taxonomy_9606_STRING.tab', root='./'
    , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))