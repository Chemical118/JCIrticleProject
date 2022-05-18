# ACIrticle Project by Julia
#### Using Random Forest in enzyme for various variables at Julia

Imported Project at [ASCIrticleProject](https://github.com/Chemical118/ASCIrticleProject) Python to Julia  
Please read [main.ipynb](https://github.com/Chemical118/JCIrticleProject/blob/master/main.ipynb) to get a example

Use only non-underbar starting function; functions start with underbar are backend only  

## Julia Module Description
`RF` : Random Forest Struct  
`RFI` : Random Forest Iteration Struct  

### Normal function
+ `nrmse` : calculate Normalized root mean square error by regr or value
+ `train_test_split` : split train data and test data

### Both `RF`, `RFI` function
+ `get_data` : using `prorf.rfclass.get_data`
+ `get_amino_loc` : using location array; `get_data(...)[2]`, return amino acid loaction array
+ `get_reg_importance` : at one condition return feature importance and variable importance
+ `iter_get_reg_importance` : iteration version of `get_reg_importance`; return mean and standard deviation
+ `view_importance` : viewer of `get_reg_importance`
+ `iter_view_importance` : viewer of `iter_get_reg_importance`
+ `rf_importance` : calculate feature importance
+ `view_alignment` : using `prorf.rfunction.view_alignment`

### Only `RFI` function
+ `get_reg_value` : RF at nfeat, ntree range return NRMSE arraay
+ `iter_get_reg_value` : iteration version of `get_reg_value`; return mean and standard deviation
+ `get_reg_value_loc` : using nfeat, ntree and NRMSE array return min loc tuple by val_mode on/off
+ `view_reg3d` : using nfeat, ntree and various NRMSE array draw 3d distribution

## PyCall Python Package Description
### `prorf.rfclass`
`RF` : Random Forest Class
+ `get_data` : return X, Y, L (location array) you must use .xls file and make index column

### `prorf.rfunction`
+ `view_alignment` : Using [Bokeh sequence aligner visualization program](https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner
), make protein sequence aligner