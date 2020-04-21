import traversalHelper as tr
from itertools import combinations
import helper as hr
import math, time
import multiprocessing
from multiprocessing import Pool, Manager
from functools import partial
import numpy as np
from statistics import mean
import json, sys, os
from collections import defaultdict

class helperFunc:
    @staticmethod
    def uvSpec_basic(node, classNodes, PPIr):
        # given a node, and a class of complementary node, the ratio of neigh node being complement over all neigh nodes
        if len(classNodes) == 0: return 0
        return len(classNodes&PPIr[node])/len(classNodes|PPIr[node])

    @staticmethod
    def uvSpec_noPad(node, classNodes, PPIr):
        if len(classNodes) == 0: return 0
        return len(classNodes&PPIr[node])/(len(classNodes|PPIr[node])-1)

    @staticmethod
    def uvSpec_linear(node, classNodes, PPIr):
        if len(classNodes) == 0: return 0
        return len(classNodes&PPIr[node])-len(PPIr[node]-classNodes)

    @staticmethod
    def xySpec_basic(node, classNodes, PPIr):
        # same as uvSpec
        if len(classNodes) == 0: return 0
        return len(classNodes&PPIr[node])/len(classNodes|PPIr[node])

    @staticmethod
    def xyContrib(node, classNodes, PPIr):
        # given nodeX (nodeY), check the ratio of degree being nodeU (nodeV)
        if len(PPIr[node]) == 0: return 0
        return len(PPIr[node]&classNodes)/len(PPIr[node])

    @staticmethod
    def uvContrib(node, parentCNodes, peerCNodes, PPIr):
        # have inherient normalization (different edge weight for node of different deg)
        if len(PPIr[node]) == 0: return 0
        return len(PPIr[node]&parentCNodes&peerCNodes)/len(PPIr[node])

    @staticmethod
    def uvContrib_padded(node, parentCNodes, peerCNodes, PPIr):
        if len(PPIr[node]) == 0: return 0
        return len(PPIr[node]&peerCNodes)/len(PPIr[node])

    @staticmethod
    def dualCN(parentNode, childNode, PPIr):
        if len(PPIr[parentNode]) == 0 and len(PPIr[childNode]) == 0: return 0
        return len(PPIr[parentNode]&PPIr[childNode])/len(PPIr[parentNode]|PPIr[childNode])

    @staticmethod
    def uvEval(parentNode, classNodes, PPIr):
        # parentNode = u or v, classNodes = V or U
        if len(PPIr[parentNode]) == 0: return 0
        return len(PPIr[parentNode]&classNodes)/len(PPIr[parentNode])

    @staticmethod
    def logging(count, lastCount, total, avgTime, startTime, frequency=1000):
        count += 1
        if count == 1: avgTime = time.time()-startTime
        else: avgTime = (avgTime*(count-1)+(time.time()-startTime))/count
        if count-lastCount > frequency:
            print("reference core's count: {}/{}. tick rate: {}. Expected sec to finish (hr): {} ({})".format(
                count, total, frequency, round(avgTime*(total-count), 2), round(avgTime*(total-count)/60/60, 2)), end="\r")
            lastCount = count
        return count, lastCount, avgTime


class ns:
    BRToRelat = tr.Helper.binary_to_relation
    toDualBR = tr.Helper.to_dual_binary_relation
    BRToNode = tr.Helper.binary_relation_to_node
    arr_pStr = tr.Helper.list_to_pathStrs
    pStr_arr = tr.Helper.pathStrs_to_list
    br_str = tr.Helper.br_to_pathStr
    L3Scoring = ["L3Normalizing", "L3uvJoin", "L3Raw"]
    L2Scoring = ["commonNeighbor"]
    CARBasedScoring = ["CRA", "CAR"]
    interStrScoring = ["interStr"]
    normFuncMapping = {"sqrt": math.sqrt, "log":math.log, "none": lambda x: x, "null": lambda x: 1}
    uvSpecMapping = {"basic": helperFunc.uvSpec_basic, "linear": helperFunc.uvSpec_linear, "noPad": helperFunc.uvSpec_noPad, "null": lambda node, classNodes, PPIr: 1}
    xyContribMapping = {"basic": helperFunc.xyContrib, "null": lambda node, classNodes, PPIr: 1}
    uvContribMapping = {"basic": helperFunc.uvContrib, "padded": helperFunc.uvContrib_padded, "null": lambda node, parentCNodes, peerCNodes, PPIr: 1}
    xySpecMapping = {"basic": helperFunc.xySpec_basic, "null": lambda node, classNodes, PPIr: 1}
    dualCNMapping = {"basic": helperFunc.dualCN, "null": lambda parentNode, childNode, PPIr: 1}
    uvJoinMapping = {"basic": True, "null": False}

def L3_normalization(PPIr, uvPair, normFunc):
    score = 0
    for uv in uvPair:
        if normFunc(len(PPIr[uv[0]])*len(PPIr[uv[1]])) == 0: continue
        score += 1/normFunc(len(PPIr[uv[0]])*len(PPIr[uv[1]]))
    return score

def get_uv(x, y, PPIr, uvJoin=False):
    candidateUs = PPIr[x]
    candidateVs = PPIr[y]
    if not uvJoin:
        candidateUs = candidateUs-candidateVs
        candidateVs = candidateVs-candidateUs
    uvPair = []
    for u in candidateUs:
        for v in candidateVs:
            if u not in PPIr[v]: continue
            uvPair.append([u,v])
    return uvPair, candidateUs, candidateVs

def deserialize_args(scoreArgs):
    normFunc, uvSpec, xySpec, xyContrib, uvContrib, dualCN, uvJoin = None, None, None, None, None, None, None
    scoreArgs = scoreArgs+["null"]*(6-len(scoreArgs))
    for i in range(len(scoreArgs)):
        if i == 0: normFunc = ns.normFuncMapping[scoreArgs[i]]
        if i == 1: uvSpec = ns.uvSpecMapping[scoreArgs[i]]
        if i == 2: xySpec = ns.xySpecMapping[scoreArgs[i]]
        if i == 3: xyContrib = ns.xyContribMapping[scoreArgs[i]]
        if i == 4: uvContrib = ns.uvContribMapping[scoreArgs[i]]
        if i == 5: dualCN = ns.dualCNMapping[scoreArgs[i]]
        if i == 6: uvJoin = ns.uvJoinMapping[scoreArgs[i]]
    return normFunc, uvSpec, xySpec, xyContrib, uvContrib, dualCN, uvJoin

def interStr_Scoring(samplePPIr, nodeX, nodeY, scoringMethod, scoreArgs):
    normFunc, uvSpec, xySpec, xyContrib, uvContrib, dualCN, uvJoin = deserialize_args(scoreArgs)
    uvPair, candidateUs, candidateVs = get_uv(nodeX, nodeY, samplePPIr, uvJoin=uvJoin)
    nodeUs, nodeVs = set([uv[0] for uv in uvPair]), ([uv[1] for uv in uvPair])
    classU, classV = set(list(nodeUs)+[nodeY]), set(list(nodeVs)+[nodeX])
    score = 0
    for [u, v] in uvPair:
        term = uvContrib(u, set([nodeX]), classV, samplePPIr)*uvContrib(v, set([nodeY]), classU, samplePPIr)
        term *= uvSpec(v, classU, samplePPIr)*uvSpec(u, classV, samplePPIr)
        term *= dualCN(nodeX, v, samplePPIr)*dualCN(nodeY, u, samplePPIr)
        score += term*1/normFunc(len(samplePPIr[u])*len(samplePPIr[v]))
    score *= xyContrib(nodeX, classU, samplePPIr)*xyContrib(nodeY, classV, samplePPIr)
    score *= xySpec(nodeX, classU, samplePPIr)*xySpec(nodeY, classV, samplePPIr)
    return score

def L3_Scoring(samplePPIr, nodeX, nodeY, scoringMethod):
    if scoringMethod == "L3Normalizing":
        uvPair, candidateUs, candidateVs = get_uv(nodeX, nodeY, samplePPIr)
        score = L3_normalization(samplePPIr, uvPair, math.sqrt)
    elif scoringMethod == "L3uvJoin":
        uvPair, candidateUs, candidateVs = get_uv(nodeX, nodeY, samplePPIr, uvJoin=True)
        score = L3_normalization(samplePPIr, uvPair, math.sqrt)
    elif scoringMethod == "L3Raw":
        uvPair, candidateUs, candidateVs = get_uv(nodeX, nodeY, samplePPIr, uvJoin=True)
        score = L3_normalization(samplePPIr, uvPair, lambda x: 1)
    return score

def L2_Scoring(samplePPIr, nodeX, nodeY, scoringMethod):
    if scoringMethod == "commonNeighbor":
        score = len(samplePPIr[nodeX]&samplePPIr[nodeY])
    return score

def CRA(samplePPIr, nodeX, nodeY):
    score, cn = 0, samplePPIr[nodeX]&samplePPIr[nodeY]
    for node in cn: score += len(samplePPIr[node]&cn)/len(samplePPIr[node])
    return score

def CAR(samplePPIr, nodeX, nodeY):
    score, cn = 0, samplePPIr[nodeX]&samplePPIr[nodeY]
    for node in cn: score += len(samplePPIr[node]&cn)/2
    return len(cn)*score

def CARBased_Scoring(samplePPIr, nodeX, nodeY, scoringMethod):
    if scoringMethod == 'CRA':
        score = CRA(samplePPIr, nodeX, nodeY)
    if scoringMethod == 'CAR':
        score = CAR(samplePPIr, nodeX, nodeY)
    return score

def _PPILinkPred(nodePairs, samplePPIr, scoringMethod, scoreArgs=[], logging=False):
    scores, predictedPPIbrs = [], []
    count, lastCount, total, avgTime = 0, 0, len(nodePairs), 0
    for nodePair in nodePairs:
        startTime = time.time()
        if nodePair[1] in samplePPIr[nodePair[0]]: continue
        nodeX, nodeY = nodePair[0], nodePair[1]
        if scoringMethod in ns.L3Scoring:
            score = L3_Scoring(samplePPIr, nodeX, nodeY, scoringMethod)
        elif scoringMethod in ns.L2Scoring:
            score = L2_Scoring(samplePPIr, nodeX, nodeY, scoringMethod)
        elif scoringMethod in ns.interStrScoring:
            score = interStr_Scoring(samplePPIr, nodeX, nodeY, scoringMethod, scoreArgs)
        elif scoringMethod in ns.CARBasedScoring:
            score = CARBased_Scoring(samplePPIr, nodeX, nodeY, scoringMethod)
        scores.append(score)
        if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)
        predictedPPIbrs.append(nodePair)
    return scores, predictedPPIbrs

def _multiCore_handler(args, iterable):
    (nodePairs, splitStartIndex, splitEndIndex, samplePPIr, scoringMethod, scoreArgs, logging, PPIresQ) = args
    nodePairs = nodePairs[splitStartIndex[iterable]:splitEndIndex[iterable]]
    logging = logging[iterable]
    scores, predictedPPIbrs = _PPILinkPred(nodePairs, samplePPIr, scoringMethod, scoreArgs, logging)
    PPIresQ.put([predictedPPIbrs, scores])
    return

def multiCore_PPILinkPred(samplePPIbr, scoringMethod, scoreArgsDict, coreNo, topNo=None, logging=False, nodePairs=None):
    # @param scoreArgs: dict, assign the normalization functions (normFunc, uvSpec, xySpec, uvContrib, xyContrib, dualCN)
    normOrder = ['normFunc', 'uvSpec', 'xySpec', 'uvContrib', 'xyContrib', 'dualCN', 'uvJoin']
    scoreArgs = ['null' if normTag not in scoreArgsDict else scoreArgsDict[normTag] for normTag in normOrder]

    samplePPIr = ns.BRToRelat(ns.toDualBR(samplePPIbr), rSet=True)
    sampleNodes = ns.BRToNode(samplePPIbr)
    if nodePairs is None: nodePairs = list(combinations(sampleNodes, 2))

    splitStartIndex = [i*math.floor(len(nodePairs)/coreNo) for i in range(0, coreNo)] # both splitting is correct
    splitEndIndex = [(i+1)*math.floor(len(nodePairs)/coreNo) if i != coreNo-1 else len(nodePairs) for i in range(0, coreNo)]
    mgr = Manager()
    PPIresQ = mgr.Queue()
    if logging: logging = [True if i == 0 else False for i in range(coreNo)]
    else: logging = [False for i in range(coreNo)]
    args = (nodePairs, splitStartIndex, splitEndIndex, samplePPIr, scoringMethod, scoreArgs, logging, PPIresQ)
    func = partial(_multiCore_handler, args)
    with Pool(coreNo) as p:
        p.map(func, [i for i in range(coreNo)])
    if logging: print("\n")
    mergedScores, mergedPPIbrs = [], []
    PPIresL = [PPIresQ.get() for i in range(coreNo)]
    for [predictedPPIbr, scores] in PPIresL:
        mergedScores += scores
        mergedPPIbrs += predictedPPIbr

    sortedPPIbrs, sortedScores = hr.sort_key_val(mergedPPIbrs, mergedScores)
    if topNo is None: topNo = len(sortedPPIbrs)
    topPredPPIbrs = sortedPPIbrs[0:topNo]
    topScores = sortedScores[0:topNo]
    return topPredPPIbrs, topScores

def get_prec(fullPPIbr, topPredPPIbrs):
    # all input supposed to be monoBR
    prec = len(set(ns.arr_pStr(topPredPPIbrs))&set(ns.arr_pStr(ns.toDualBR(fullPPIbr))))/len(topPredPPIbrs)
    return prec
    
def get_sliding_prec(fullPPIbr, topPredPPIbrs, loopRange, logging):
    minI, maxI = loopRange[0], loopRange[1]
    fullPPIbr = set(ns.arr_pStr(ns.toDualBR(fullPPIbr)))
    precTop, precBot = len(set(ns.arr_pStr(topPredPPIbrs[0:minI]))&fullPPIbr), minI
    precTops = []
    precs = [precTop/precBot]
    count, lastCount, total, avgTime, startTime = minI+1, 0, maxI, 0, time.time()
    for i in range(minI+1, maxI):
        if ns.br_str(topPredPPIbrs[i]) in fullPPIbr or ns.br_str(topPredPPIbrs[i][::-1]) in fullPPIbr: precTop += 1
        precTops.append(precTop)
        if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)
    precTops = np.asarray(precTops)
    precBots = np.asarray([i for i in range(minI+1, maxI)])
    precs += list(np.divide(precTops, precBots))
    return precs

def get_rec(fullPPIbr, samplePPIbrs, topPredPPIbrs):
    relevant = ns.pStr_arr(set(ns.arr_pStr(fullPPIbr))-set(ns.arr_pStr(ns.toDualBR(samplePPIbrs))))
    rec = len(set(ns.arr_pStr(ns.toDualBR(topPredPPIbrs)))&set(ns.arr_pStr(fullPPIbr)))/len(relevant)
    return rec

def get_sliding_rec(fullPPIbr, samplePPIbrs, topPredPPIbrs, loopRange, logging):
    minI, maxI = loopRange[0], loopRange[1]
    relevant = set(ns.arr_pStr(fullPPIbr))-set(ns.arr_pStr(ns.toDualBR(samplePPIbrs)))
    fullPPIbr = set(ns.arr_pStr(ns.toDualBR(fullPPIbr)))
    recTop, recBot = len(set(ns.arr_pStr(topPredPPIbrs[0:minI]))&fullPPIbr), len(relevant)
    recTops = []
    recs = [recTop/recBot]
    count, lastCount, total, avgTime, startTime = minI+1, 0, maxI, 0, time.time()
    for i in range(minI+1, maxI):
        if ns.br_str(topPredPPIbrs[i]) in fullPPIbr or ns.br_str(topPredPPIbrs[i][::-1]) in fullPPIbr: recTop += 1
        recTops.append(recTop)
        if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)
    recTops = np.asarray(recTops)
    recBots = np.asarray([recBot for i in range(minI+1, maxI)])
    recs += list(np.divide(recTops, recBots))
    return recs

def precRecMap_handler(args, iterable):
    (perTagNo, tags, predPPIbr, samplePPIbr, fullPPIbr, logging, resQ) = args
    logging = logging[iterable]
    partialPrecRecMap = {}
    for tagNo in perTagNo[iterable]:
        loopRange = (1, len(predPPIbr[tagNo]))
        partialPrecRecMap[tags[tagNo]] = {}
        partialPrecRecMap[tags[tagNo]]["prec"] = get_sliding_prec(fullPPIbr[tagNo], predPPIbr[tagNo], loopRange, logging)
        partialPrecRecMap[tags[tagNo]]["rec"] = get_sliding_rec(fullPPIbr[tagNo], samplePPIbr[tagNo], predPPIbr[tagNo], loopRange, logging)
    resQ.put(partialPrecRecMap)
    return

def precRecMap_multiCore(tags, predPPIbr, samplePPIbr, fullPPIbr, coreNo, logging=False):
    tagNo, perTagNo = [i for i in range(len(tags))], [[] for i in range(coreNo)]
    for i in range(math.ceil(len(tags)/coreNo)):
        for j in range(coreNo*i, min(coreNo*(i+1), coreNo*i+(len(tags)-coreNo*i))):
            perTagNo[j-coreNo*i].append(tagNo[j])
    mgr = Manager()
    resQ = mgr.Queue()
    if logging: logging = [True if i == 0 else False for i in range(coreNo)]
    else: logging = [False for i in range(coreNo)]
    args = (perTagNo, tags, predPPIbr, samplePPIbr, fullPPIbr, logging, resQ)
    func = partial(precRecMap_handler, args)
    with Pool(coreNo) as p:
        p.map(func, [i for i in range(coreNo)])
    precRecMap = {}
    for i in range(coreNo): precRecMap.update(resQ.get())
    return precRecMap

if __name__ == "__main__":
    pass