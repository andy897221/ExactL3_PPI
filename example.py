import sys
sys.path.append("./src")
import src.onPPILinkPred as ppiLPred
from itertools import combinations
import src.helper as helper

# consider such a PPI network
#  A
# / \
#B  E
#|  |
#C  F
#\ /
# D
# the binary relation of such ppi networks can be expressed as such a dictionary
PPIr = {"A": {"B", "E"}, "B": {"A", "C"}, "C": {"B", "D"}
, "E": {"A", "F"}, "F": {"E", "D"}, "D": {"C", "F"}}
# or such a list
nodes = ["A", "B", "C", "D", "E", "F"]
# we then generate all possible edges (the algorithm will filter non-candidate edges)
nodePairs = list(combinations(nodes, 2))

# to generate all candidate edges with scores for the original L3 algorithm
L3PPIs, L3Scores = ppiLPred._PPILinkPred(nodePairs, PPIr, scoringMethod="L3uvJoin")
L3sortedPPIs, L3sortedScores = helper.sort_key_val(L3PPIs, L3Scores)
print(L3sortedPPIs)
print(L3sortedScores)
print("")

# For ExactL3_1, generate as such.
# we defined the simple penalization index as xyContrib in ppiLPred.helperFunc, and jaccard coefficient as dualCN
# so, we pass this as arguments via our custom scoreArgs parameter
# the scoreArg is arranged as a list, ['normFunc'=?, 'uvSpec'=?, 'xySpec'=?, 'uvContrib'=?, 'xyContrib'=?, 'dualCN'=?, 'uvJoin'=?]
# for our case, it would be ['null', 'null', 'null', 'null', 'basic', 'basic', 'basic']
ExactL3_1_PPIs, ExactL3_1_Scores = ppiLPred._PPILinkPred(nodePairs, PPIr, scoringMethod="interStr", 
    scoreArgs=['null', 'null', 'null', 'null', 'basic', 'basic', 'basic'])
ExactL3_1_sortedPPIs, ExactL3_1_sortedScores = helper.sort_key_val(ExactL3_1_PPIs, ExactL3_1_Scores)
print(ExactL3_1_sortedPPIs)
print(ExactL3_1_sortedScores)
print("")

# For ExactL3_2, generate as such.
# ['null', 'basic', 'basic', 'null', 'basic', 'null', 'basic']
# note that xyContrib = xySpec in this case, since U is a subset of n(X), so U union n(X) = n(X)
ExactL3_2_PPIs, ExactL3_2_Scores = ppiLPred._PPILinkPred(nodePairs, PPIr, scoringMethod="interStr", 
    scoreArgs=['null', 'basic', 'basic', 'null', 'basic', 'null', 'basic'])
ExactL3_2_sortedPPIs, ExactL3_2_sortedScores = helper.sort_key_val(ExactL3_2_PPIs, ExactL3_2_Scores)
print(ExactL3_2_sortedPPIs)
print(ExactL3_2_sortedScores)
print("")


# example for accessing the database
# available dataset: bioGRID, STRING, MINT, IntAct
import bioGRID as bg
bg_GI, bg_PPI = bg.parse_bioGRID(root="./src")
print(bg_PPI.head())
print(bg_GI.head())

# for other datasets:
# import STRING, MINT, IntAct
# STRING.parse_STRING(root="./src")
# MINT.parse_MINT(root="./src")
# IntAct.parse_IntAct(spokeModel=True, root="./src")