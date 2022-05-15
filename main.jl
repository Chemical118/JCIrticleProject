include("prorf.jl")
import .prorf as pf

using DecisionTree

R = pf.RF("Data/rgpdata.fasta", "Data/rdata.xls")
X, Y, L = pf.get_data(R, 12, 'B')

pf.get_reg_importance(R, X, Y, L, 5, 400)