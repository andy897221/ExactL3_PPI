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
    # only direct interaction & physical association & association are PPIs
    df = pd.read_csv(ppiFile, sep='\t', header=None)

    ppi_df = df[df[11].isin([i for i in df[11] if 'direct interaction' in i or 'physical association' in i or 'association' in i])].copy()
    inv_ppis = np.transpose(np.asarray([ppi_df[0], ppi_df[1]]))
    ppis = []
    for i in inv_ppis:
        if i[0] == "-" or i[1] == "-": continue
        ppi_h = "-".join(i[0].split(":")[1].split("-")[:-1]) if len(i[0].split(":")[1].split("-")) > 1 else i[0].split(":")[1].split("-")[0]
        ppi_t = "-".join(i[1].split(":")[1].split("-")[:-1]) if len(i[1].split(":")[1].split("-")) > 1 else i[1].split(":")[1].split("-")[0]
        ppis.append([ppi_h, ppi_t])

    mappedPPIs = []
    for ppi in ppis:
        if ppi[0] not in uniProtMapping or ppi[1] not in uniProtMapping: continue
        mappedPPIs.append([uniProtMapping[ppi[0]], uniProtMapping[ppi[1]]])
    mappedPPIs = np.transpose(np.asarray(mappedPPIs))
    ppi_df = pd.DataFrame({'nodeA': mappedPPIs[0], 'nodeB': mappedPPIs[1]})


    ggi_df = df[~df[11].isin([i for i in df[11] if 'direct interaction' in i or 'physical association' in i or 'association' in i])].copy()
    inv_ppis = np.transpose(np.asarray([ggi_df[0], ggi_df[1]]))
    ppis = []
    for i in inv_ppis:
        if i[0] == "-" or i[1] == "-": continue
        ppi_h = "-".join(i[0].split(":")[1].split("-")[:-1]) if len(i[0].split(":")[1].split("-")) > 1 else i[0].split(":")[1].split("-")[0]
        ppi_t = "-".join(i[1].split(":")[1].split("-")[:-1]) if len(i[1].split(":")[1].split("-")) > 1 else i[1].split(":")[1].split("-")[0]
        ppis.append([ppi_h, ppi_t])

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

    ggi_df, ppi_df = parse_MINT(ppiFile='./data/MINT/species human', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
                wFile_GGI='./data/parsed/MINT_homo_GGI.pkl', wFile_PPI='./data/parsed/MINT_homo_PPI.pkl', root='./')
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))