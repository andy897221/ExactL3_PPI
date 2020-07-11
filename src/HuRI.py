import pandas as pd
import numpy as np
import os

def parse_uniProt_map(uniProtMap):
    df = pd.read_csv(uniProtMap, sep='\t')
    df.dropna(inplace=True)
    uniProtMapping = dict(zip(list(df['Entry']), list(df['Gene names  (primary )'])))
    return uniProtMapping

def parse_HuRI(ppiFile='./data/atlas/HuRI.psi', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
                wFile_PPI='./data/parsed/HuRI_PPI.pkl', root='./'):
    ppiFile, uniProtMap, wFile_PPI = root+ppiFile, root+uniProtMap, root+wFile_PPI
    if os.path.exists(wFile_PPI): return pd.read_pickle(wFile_PPI)

    uniProtMapping = parse_uniProt_map(uniProtMap)
    # only direct interaction & physical association & association are PPIs
    ppi_df = pd.read_csv(ppiFile, sep='\t', header=None)
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

    ppi_df.to_pickle(wFile_PPI)
    return ppi_df


if __name__ == "__main__":
    ppi_df = parse_HuRI()
    print(ppi_df.head())
    print(len(ppi_df.index))