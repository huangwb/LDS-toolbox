function [alpha,D,qX]=sparse_coding(data,dict,SC_opts)
% Sparse coding on the space of linear dynamical systems
%
% INUTS
% data            - {} LDS data
% dict            - {} LDS dictionary
% beta            - [] The weight determines the trade-off between the mean component and covariance component. Setting 1 means no covariance involved, while 0 means no mean involved. see demo.m
% SC_opts      - .  SC set-up parameters, see demo.m
%
% OUTPUTS
% alpha        - [] Returns the sparse codes
% D            - [] Ditionary kernel matrix
% qX           - [] Dictioanry-Data kernel matrix
% implemented by Wenbing Huanng, 2015-3-23

Solver_Flag = 1;

nAtoms = size(dict,2);
nPoints = size(data,2);

computeKernel = SC_opts.kernel;

if isfield(SC_opts, 'K_D')
    K_D = SC_opts.K_D;
else
    K_D = computeKernel(dict);
end
K_XD =computeKernel(dict,data);

[KD_U,KD_D,~] = svd(K_D);
D = diag(sqrt(diag(KD_D)))*KD_U';
D_Inv = KD_U*diag(1./sqrt(diag(KD_D)));
qX = D_Inv'*K_XD;

switch Solver_Flag
    case 1
        alpha = full(mexLasso(qX,D,struct('mode',2,'lambda',SC_opts.lambda,'lambda2',0)));
    otherwise 
        alpha = zeros(nAtoms,nPoints);       
        for t1 = 1:nPoints 
            cvx_begin quiet;
            variable alpha_cvx(nAtoms,1);
            minimize( norm(qX(:,t1) - D * alpha_cvx) +  SC_opts.lambda*norm(alpha_cvx,1));
            cvx_end;
            alpha(:,t1) = double(alpha_cvx);
        end  
end

end