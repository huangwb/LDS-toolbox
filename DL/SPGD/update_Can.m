function newDict = update_Can(data,dict,alpha,SPGD_opts)
%UPDATE_CAN Update the dictionary under canonical representation
%
% INUTS
% data            - {} LDS data
% dict            - {} LDS dictionary
% alpha           - [] Sparse codes
%
% OUTPUTS
% newDict        - {} Updated dictionary
% implemented by Wenbing Huanng, 2016-12-12

MIN_THRESH = 1e-6;
nAtoms = length(dict);
nData = length(data);
n = size(dict{1}.S,1);

for r=1:nAtoms
    
    alpha_nonzero = find(alpha(r,:),1);
    if (isempty(alpha_nonzero))
%         fprintf('a useless atom identified!\n');
        continue;
    end
    delta_L_U = 0;
    delta_L_S = 0;
    delta_L_V = 0;
    
    A_r = diag(dict{r}.S)*dict{r}.V; 
    S_rr = dlyap(A_r',eye(n));
    [tmpU, tmpS, tmpV] = svd(S_rr);
    iS_rr = tmpV*diag(1./diag(tmpS))*tmpU';
    
    % compute gradients over dictionary
    for id_atom=1:nAtoms
        tmpAlpha_rj = alpha(r,:)*alpha(id_atom,:)';
        if id_atom == r || abs(tmpAlpha_rj) < MIN_THRESH
            continue;
        end
        A_j = diag(dict{id_atom}.S)*dict{id_atom}.V;
        
        S_rj = dlyap(A_r',A_j,dict{r}.U'*dict{id_atom}.U);
        S_jj = dlyap(A_j',eye(n));
        
        [tmpU, tmpS, tmpV] = svd(S_jj);
        iS_jj = tmpV*diag(1./diag(tmpS))*tmpU';        
        
        G_rr =  -iS_rr*S_rj*iS_jj*S_rj'*iS_rr;
        G_rj = iS_rr*S_rj*iS_jj;
        
        R_rr = dlyap(A_r,G_rr);
        R_rj = dlyap(A_r, A_j',G_rj);
        
        delta_L_A = 2*tmpAlpha_rj*(S_rr*A_r*R_rr + S_rj*A_j*R_rj');
        delta_L_S = delta_L_S + diag(delta_L_A*dict{r}.V');
        delta_L_V = delta_L_V + diag(dict{r}.S)*delta_L_A; 
        
        delta_L_U = delta_L_U + 2*tmpAlpha_rj*(dict{r}.U*R_rr + dict{id_atom}.U*R_rj');
        
%         if ~isreal(delta_L_S) || (r==17) 
%             stop = 1;
%         end
        
    end
    
    % compute gradients over data
    for id_data=1:nData
        if alpha(r,id_data)< MIN_THRESH
            continue;
        end
        A_i = diag(data{id_data}.S)*data{id_data}.V;
        
        S_ri = dlyap(A_r',A_i,dict{r}.U'*data{id_data}.U);
        S_ii = dlyap(A_i',eye(n));
        
        [tmpU, tmpS, tmpV] = svd(S_ii);
        iS_ii = tmpV*diag(1./diag(tmpS))*tmpU';        
        
        G_rr =  -iS_rr*S_ri*iS_ii*S_ri'*iS_rr;
        G_ri = iS_rr*S_ri*iS_ii;
        
        R_rr = dlyap(A_r,G_rr);
        R_ri = dlyap(A_r, A_i',G_ri);
        
        delta_L_A = -2*alpha(r,id_data)*(S_rr*A_r*R_rr + S_ri*A_i*R_ri');
        delta_L_S = delta_L_S + diag(delta_L_A*dict{r}.V');
        delta_L_V = delta_L_V + diag(dict{r}.S)*delta_L_A;
        
        delta_L_U = delta_L_U - 2*alpha(r,id_data)*(dict{r}.U*R_rr + data{id_data}.U*R_ri');
        
    end
    
    % normalize the gradients
    delta_L_U = 1/nData * delta_L_U;
    delta_L_S = 1/nData * delta_L_S;
    delta_L_V = 1/nData * delta_L_V;
    
    % update S
    delta_S = SPGD_opts.learningRate_A;
    dict{r}.S = dict{r}.S - delta_S * delta_L_S;
    dict{r}.S = 1./max(abs(dict{r}.S),1) .* dict{r}.S;
%     a = SPGD_opts.SN_a;
%     b = SPGD_opts.SN_b;
%     rho_old = 1/a*(log(b+dict{r}.S)-log(b-dict{r}.S));
%     rho_new = rho_old -2*delta_S*a*b*delta_L_S.*sigm(a*rho_old).*(1-sigm(a*rho_old));
%     dict{r}.S = real(b*(2*sigm(a*rho_new)-1)); % in case that rho_old and rho_new are complex when previous dict{r}.S=1
    
    % update V
    delta_V = SPGD_opts.learningRate_A;
    dict{r}.V = project_to_manifold(dict{r}.V, delta_L_V, delta_V);
    
    % update U
    delta_t = SPGD_opts.learningRate_C;
    dict{r}.U = project_to_manifold(dict{r}.U, delta_L_U, delta_t);
    
    % update iL
    A_r = diag(dict{r}.S)*dict{r}.V; 
    S_rr = dlyap(A_r',eye(n));
    [Udot,Sdot,~] = svd(S_rr);
    dict{r}.iL = diag(1./sqrt(diag(Sdot)))*Udot';
    
end


newDict = dict;

end

