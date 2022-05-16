include("prorf.jl")
import .prorf as pf

R = pf.RF("Data/rgpdata.fasta", "Data/rdata.xls")
X, Y, L = pf.get_data(R, 9, 'B')

M, Z = pf.get_reg_importance(R, X, Y, L, 5, 400, val_mode=true)