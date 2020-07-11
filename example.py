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

# For ExactL3, generate as such.
# we defined the simple penalization index as xyContrib in ppiLPred.helperFunc, and jaccard coefficient as dualCN
# so, we pass this as arguments via our custom scoreArgs parameter
# the scoreArg is arranged as a list, ['normFunc'=?, 'uvSpec'=?, 'xySpec'=?, 'uvContrib'=?, 'xyContrib'=?, 'dualCN'=?, 'uvJoin'=?]
# for our case, it would be ['null', 'null', 'null', 'null', 'basic', 'basic', 'basic']
ExactL3_PPIs, ExactL3_Scores = ppiLPred._PPILinkPred(nodePairs, PPIr, scoringMethod="interStr", 
    scoreArgs=['null', 'null', 'null', 'null', 'basic', 'basic', 'basic'])
ExactL3_sortedPPIs, ExactL3_sortedScores = helper.sort_key_val(ExactL3_PPIs, ExactL3_Scores)
print(ExactL3_sortedPPIs)
print(ExactL3_sortedScores)
print("")


# example for accessing the database
# available dataset: bioGRID, STRING, MINT, IntAct, HuRI
import bioGRID as bg
bg_GI, bg_PPI = bg.parse_bioGRID(root="./src")
print(bg_PPI.head())
print(bg_GI.head())

# for other datasets:
# import STRING, MINT, IntAct
# STRING.parse_STRING(root="./src")
# MINT.parse_MINT(root="./src")
# IntAct.parse_IntAct(spokeModel=True, root="./src")
# HuRI.parse_HuRI(root="./src")
# by default they all returns Yeast dataset, please see the notebook and other provided data generating script to return the Human dataset, or just ask me