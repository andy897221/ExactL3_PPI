import pandas as pd
import numpy as np
import os

def parse_uniProt_map(uniProtMap):
    df = pd.read_csv(uniProtMap, sep='\t')
    df.dropna(inplace=True)
    uniProtMapping = dict(zip(list(df['Entry']), list(df['Gene names  (primary )'])))
    return uniProtMapping

def parse_MINT(ppiFile='./data/MINT/species yeast', uniProtMap="./data/UniProt/uniprot-taxonomy_559292.tab", 
                wFile_GGI='./data/parsed/MINT_GGI.pkl', wFile_PPI='./data/parsed/MINT_PPI.pkl', root='./'):
    ppiFile, uniProtMap, wFile_GGI, wFile_PPI = root+ppiFile, root+uniProtMap, root+wFile_GGI, root+wFile_PPI
    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)

    uniProtMapping = parse_uniProt_map(uniProtMap)
    # only direct interaction & physical association & colocalization & association are PPIs
    df = pd.read_csv(ppiFile, sep='\t', header=None)
    
    ppi_df = df[df[11].isin([i for i in df[11] if 'direct interaction' in i or 'physical association' in i or 'colocalization' in i or 'association' in i])].copy()
    ppis = np.transpose(np.asarray([[i.split(":")[1] for i in ppi_df[0]], [i.split(":")[1] for i in ppi_df[1]]]))
    mappedPPIs = []
    for ppi in ppis:
        if ppi[0] not in uniProtMapping or ppi[1] not in uniProtMapping: continue
        mappedPPIs.append([uniProtMapping[ppi[0]], uniProtMapping[ppi[1]]])
    mappedPPIs = np.transpose(np.asarray(mappedPPIs))
    ppi_df = pd.DataFrame({'nodeA': mappedPPIs[0], 'nodeB': mappedPPIs[1]})

    ggi_df = df[~df[11].isin([i for i in df[11] if 'direct interaction' in i or 'physical association' in i or 'colocalization' in i or 'association' in i])].copy()
    ppis = np.transpose(np.asarray([[i.split(":")[1] for i in ggi_df[0]], [i.split(":")[1] for i in ggi_df[1]]]))
    mappedPPIs = []
    for ppi in ppis:
        if ppi[0] not in uniProtMapping or ppi[1] not in uniProtMapping: continue
        mappedPPIs.append([uniProtMapping[ppi[0]], uniProtMapping[ppi[1]]])
    mappedPPIs = np.transpose(np.asarray(mappedPPIs))
    ggi_df = pd.DataFrame({'nodeA': mappedPPIs[0], 'nodeB': mappedPPIs[1]})

    ppi_df.to_pickle(wFile_PPI)
    ggi_df.to_pickle(wFile_GGI)
    return ggi_df, ppi_df


if __name__ == "__main__":
    ggi_df, ppi_df = parse_MINT()
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))