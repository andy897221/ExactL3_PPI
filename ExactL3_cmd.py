if __name__ == "__main__":
    import sys
    sys.path.append("./src")
    import src.onPPILinkPred as ppiLPred
    from itertools import combinations
    import src.helper as helper
    import pandas as pd
    import numpy as np

    ppiFile, savePath, method, coreNo = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

    ppi_df = pd.read_csv(ppiFile, sep="\t", header=None)
    ppi_br = ppi_df[[0,1]].values.tolist()

    scoreArgs = {}
    if method == "ExactL3":
        scoreArgs = {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}
        algor = "interStr"
    elif method == "L3":
        algor = "L3uvJoin"
    elif method == "CN":
        algor = "commonNeighbor"
    elif method == "CH2_L3":
        algor = "CH2_L3"
    elif method == "Sim":
        algor = "Sim"
    elif method == "CRA":
        algor = "CRA"
    else:
        print("link predictor '{}' is invalid.".format(method))
        exit()

    sortedPPIs, sortedScores = ppiLPred.multiCore_PPILinkPred(ppi_br, algor, scoreArgs, coreNo)
    sortedPPIs = np.asarray(sortedPPIs).transpose()
    save_df = pd.DataFrame({"proteinA": sortedPPIs[0], "proteinB": sortedPPIs[1], "score": sortedScores})
    save_df.to_csv(savePath, sep="\t", index=False)