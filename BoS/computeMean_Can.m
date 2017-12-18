function center = computeMean_Can(data, opts, varargin)
%COMPUTMEAN_MDS Summary of this function goes here
%   Detailed explanation goes here

n = length(data{1}.S);
nData = length(data);
computeKernel = opts.kernel;

% initialization
if nargin < 3
	% compute kernel matrix
	distM = computeKernel(data)/n;
	distM = distM - diag(diag(distM)) + diag(ones(nData,1));
	[~,minIdx] = max(sum(distM));
	center = data{minIdx};
else
	center = varargin{1};
end


for iter=1:opts.nIter
    
    delta_L_U = 0;
    delta_L_S = 0;
    delta_L_V = 0;
    
    A_r = diag(center.S)*center.V;
    S_rr = dlyap(A_r',eye(n));
    
    [tmpU, tmpS, tmpV] = svd(S_rr);
    iS_rr = tmpV*diag(1./diag(tmpS))*tmpU';
    
    % compute gradients over data
    for id_data=1:nData
        
        A_i = diag(data{id_data}.S)*data{id_data}.V;
        
        S_ri = dlyap(A_r',A_i,center.U'*data{id_data}.U);
        S_ii = dlyap(A_i',eye(n));
        
        [tmpU, tmpS, tmpV] = svd(S_ii);
        iS_ii = tmpV*diag(1./diag(tmpS))*tmpU';
        
        G_rr =  -iS_rr*S_ri*iS_ii*S_ri'*iS_rr;
        G_ri = iS_rr*S_ri*iS_ii;
        
        R_rr = dlyap(A_r,G_rr);
        R_ri = dlyap(A_r, A_i',G_ri);
        
        delta_L_A = -2*(S_rr*A_r*R_rr + S_ri*A_i*R_ri');
        delta_L_S = delta_L_S + diag(delta_L_A*center.V');
        delta_L_V = delta_L_V + diag(center.S)*delta_L_A;
        
        delta_L_U = delta_L_U - 2*(center.U*R_rr + data{id_data}.U*R_ri');
       
    end
    
    % normalize the gradients
    delta_L_U = 1/nData * delta_L_U;
    delta_L_S = 1/nData * delta_L_S;
    delta_L_V = 1/nData * delta_L_V;
    
    % update S
    delta_S = opts.learningRate_A;
    center.S = center.S - delta_S * delta_L_S;
    center.S = 0.999/max(max(abs(center.S)),0.999) * center.S;
%     delta_S = opts.learningRate_A;
%     a = 2.5;
%     b = 1;
%     rho_old = 1/a*(log(b+center.S)-log(b-center.S));
%     rho_new = rho_old -2*delta_S*a*b*delta_L_S.*sigm(a*rho_old).*(1-sigm(a*rho_old));
%     center.S = real(b*(2*sigm(a*rho_new)-1)); % in case that rho_old and rho_new are complex when previously center.S=1
    
    % update V
    delta_V = opts.learningRate_A;
    center.V = project_to_manifold(center.V, delta_L_V, delta_V);
    
    % update U
    delta_t = opts.learningRate_C;
    center.U = project_to_manifold(center.U, delta_L_U, delta_t);

    % update iL
    A_r = diag(center.S)*center.V; 
    S_rr = dlyap(A_r',eye(n));
    [Udot,Sdot,~] = svd(S_rr);
    center.iL = diag(1./sqrt(diag(Sdot)))*Udot';

end


end

