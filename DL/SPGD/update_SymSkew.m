function newDict = update_SymSkew(data,dict,alpha,SPGD_opts)
%UPDATE_SYMSKEW Update the two-fold dictionary via SGD
%
% INUTS
% data            - {} LDS data
% dict            - {} LDS dictionary
% alpha           - [] Sparse codes
%
% OUTPUTS
% newDict        - {} Updated dictionary
% implemented by Wenbing Huanng, 2016-6-12

newDict = update_Sym(data,dict,alpha,SPGD_opts);
newDict = update_Skew(data,newDict,alpha,SPGD_opts);

end

