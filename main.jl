include("prorf.jl")
import .prorf as pf

R = pf.RFI("Data/rgpdata.fasta", "Data/rdata.xls", 2:1:4, 50:100:500)
X, Y, L = pf.get_data(R, 9, 'C')

pf.get_reg_value(R, X, Y)