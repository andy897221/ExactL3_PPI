import os
import pandas as pd
import traversalHelper as tr
import numpy as np
from collections import defaultdict
from statistics import mean
pd.options.mode.chained_assignment = None

# data standard: nodeA, nodeB, type, score

def parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.176.tab2.txt', root='./'
    , wFile_GGI='./data/parsed/BioGRID_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_PPI.pkl'):
    wFile_GGI, wFile_PPI, filename = root+wFile_GGI, root+wFile_PPI, root+filename

    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)

    bioGRID_df = pd.read_csv(filename, sep='\t')
    colMap = {"Official Symbol Interactor A": "nodeA",
     "Official Symbol Interactor B": 'nodeB',
      "Experimental System Type": "type", "Score": "score"}
    
    ggi_df = bioGRID_df[bioGRID_df["Experimental System Type"] == 'genetic']
    ppi_df = bioGRID_df[bioGRID_df["Experimental System Type"] == 'physical']

    for df in [ggi_df, ppi_df]:
        for col in df.columns:
            if col not in colMap: df.drop(col, axis=1, inplace=True)
        df.rename(colMap, axis=1, inplace=True)
        df.reset_index(inplace=True)

    for df in [ggi_df, ppi_df]:
        sortedBR = tr.Helper.list_to_pathStrs([sorted(br) for br in np.asarray(df[['nodeA', 'nodeB']])])
        scoreMap, typeMap = defaultdict(list), {}
        for i in range(len(sortedBR)): # remove directional interaction
            if df['nodeA'][i] == df['nodeB'][i]: continue # remove self-interaction
            scoreMap[sortedBR[i]].append(df['score'][i])
            typeMap[sortedBR[i]] = df['type'][i]
        for k, arr in scoreMap.items():
            if '-' in arr and len(arr) == 2: arr.remove('-')
            elif len(arr) == 2: arr = [mean([float(i) for i in arr])]
            arr = arr[0]
        ppi = list(np.transpose(np.asarray(tr.Helper.pathStrs_to_list(list(scoreMap.keys())))))
        df = pd.DataFrame({'nodeA': ppi[0], 'nodeB': ppi[1],
         'type': list(typeMap.values()), 'score': list(scoreMap.values())})

    ggi_df.to_pickle(wFile_GGI)
    ppi_df.to_pickle(wFile_PPI)
    return ggi_df, ppi_df

def parse_geneName_map(filename="./data/BioGRID/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.176.tab2.txt", root="./"):
    #  systematic name => gene name
    bioGRID_df = pd.read_csv(root+filename, sep="\t")
    sysName = list(bioGRID_df['Systematic Name Interactor A'])+list(bioGRID_df['Systematic Name Interactor B'])
    geneName = list(bioGRID_df['Official Symbol Interactor A'])+list(bioGRID_df['Official Symbol Interactor B'])
    geneMap = dict(zip(geneName, sysName))
    reverseGeneMap = dict(zip(sysName, geneName))
    return geneMap, reverseGeneMap

if __name__ == "__main__":
    parse_bioGRID()