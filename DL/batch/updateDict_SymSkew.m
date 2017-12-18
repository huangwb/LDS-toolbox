function newDict = updateDict_SymSkew(data,dict,alpha,LDS_opts)

newDict = updateDict_Sym(data,dict,alpha,LDS_opts);
newDict = updateDict_Skew(data,newDict,alpha,LDS_opts);


end

