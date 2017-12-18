function [alpha,D,qX]=local_sparse_coding(data,dict,SC_opts)
% local sparse coding on the space of linear dynamical systems
% INUTS
% data            - {} LDS data
% dict            - {} LDS dictionary
% SC_opts      - .  SC set-up parameters, see ./params/param_default
%
% OUTPUTS
% alpha        - [] Returns the sparse codes
% D            - [] Ditionary kernel matrix
% qX           - [] Dictioanry-Data kernel matrix
% implemented by Wenbing Huanng, 2015-3-23

computeKernel = SC_opts.kernel;

if isfield(SC_opts, 'K_D')
    K_D = SC_opts.K_D;
else
    K_D = computeKernel(dict);
end
K_XD =computeKernel(dict,data);

if isfield(data{1},'A')
   p = size(data{1}.A,1);
elseif isfield(data{1},'A1')
    p = size(data{1}.A1,1);
elseif isfield(data{1}, 'S');
    p = size(data{1}.S,1);
elseif isfield(data{1},'O')
    p = size(data{1}.O,2);
end
% % p = size(data{1}.A,2);
% q = size(data{1}.R,2);
% 
% computeKernel = SC_opts.kernel;
% 
K_D = 1/p*K_D;
K_XD =1/p*K_XD;


[KD_U,KD_D,~] = svd(K_D);
D = diag(sqrt(diag(KD_D)))*KD_U';
D_Inv = KD_U*diag(1./sqrt(diag(KD_D)));
qX = D_Inv'*K_XD;

alpha = full(mexOMP(qX,D,struct('L',SC_opts.L)));

end