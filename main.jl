include("prorf.jl")
import .prorf as pf

using DecisionTree

R = pf.RF("Data/rgpdata.fasta", "Data/rdata.xls")
X, Y, L = pf.get_data(R, 12, 'B')

x_train, x_test, y_train, y_test = pf.train_test_split(X, Y, test_size=0.7, random_state=100)
