import onPPILinkPred as ExactL3
import traversalHelper as tr

scoreArgsDictList = {
    "ExactL3_1": {'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'},
    "ExactL3_2": {'xyContrib': 'basic', 'uvSpec': 'basic', 'xySpec': 'basic', 'uvJoin': 'basic'}
}

def toScoreArgs(scoreArgsDict):
    normOrder = ['normFunc', 'uvSpec', 'xySpec', 'uvContrib', 'xyContrib', 'dualCN', 'uvJoin']
    scoreArgs = ['null' if normTag not in scoreArgsDict else scoreArgsDict[normTag] for normTag in normOrder]
    return scoreArgs

# gen xy, then 10 pairs of uv
U = [i for i in range(0, 10)]
V = [i for i in range(10, 20)]
X, Y = 20, 21

# 0 to 10, each connect 1 path of length, connect v to x first exhaustively, then u to y
for a in range(0, 20):
    # build BR for all P_4
    br = [(U[i], V[i]) for i in range(len(U))] + [(X,u) for u in U] + [(Y, v) for v in V]
    # build n P_3
    for i in range(a):
        if i >= len(V): br.append([Y, U[i-len(V)]])
        else: br.append([X, V[i]])
    samplePPIr = tr.Helper.binary_to_relation(tr.Helper.to_dual_binary_relation(br), rSet=True)
    score1, _ = ExactL3._PPILinkPred([[str(X),str(Y)]], samplePPIr, scoringMethod="interStr", scoreArgs=toScoreArgs(scoreArgsDictList["ExactL3_1"]), logging=False)
    score2, _ = ExactL3._PPILinkPred([[str(X),str(Y)]], samplePPIr, scoringMethod="interStr", scoreArgs=toScoreArgs(scoreArgsDictList["ExactL3_2"]), logging=False)
    score3, _ = ExactL3._PPILinkPred([[str(X),str(Y)]], samplePPIr, scoringMethod="L3uvJoin", scoreArgs={}, logging=False)

    print("no. of P_3: {}, ExactL3_1 score: {}".format(a, "{:5f}".format(score1[0])))
    print("no. of P_3: {}, ExactL3_2 score: {}".format(a, "{:5f}".format(score2[0])))
    print("no. of P_3: {}, L3 score: {}".format(a, "{:5f}".format(score3[0])))
    print("")