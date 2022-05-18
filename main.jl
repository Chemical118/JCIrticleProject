include("prorf.jl")
import .prorf as pf
using BenchmarkTools

R = pf.RFI("Data/rgpdata.fasta", "Data/rdata.xls", 2:1:4, 50:100:250)
X, Y, L = pf.get_data(R, 9, 'D')

MF, SF = pf.iter_get_reg_importance(R, X, Y, L, 3, 100, 10);