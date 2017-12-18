function newDict = update_Sym(data,dict,alpha,SPGD_opts)
%SGD_SYM Update the symmetric dictionary via SGD
%
% INUTS
% data            - {} LDS data
% dict            - {} LDS dictionary
% alpha           - [] Sparse codes
%
% OUTPUTS
% newDict        - {} Updated dictionary
% implemented by Wenbing Huanng, 2016-6-12

MIN_THRESH = 1e-6;
nAtoms = length(dict);
nData = length(data);
n = size(dict{1}.A1,1);
%------ update A1 and C1 -------
for r=1:nAtoms
    
    alpha_nonzero = find(alpha(r,:),1);
    if (isempty(alpha_nonzero))
        fprintf('a useless atom identified!\n');
        continue;
    end
    delta_L_C1 = 0;
    delta_L_A1 = 0;
    
    % compute gradients over dictionary
    for id_atom=1:nAtoms
        tmpAlpha_rj = alpha(r,:)*alpha(id_atom,:)';
        if id_atom == r || abs(tmpAlpha_rj) < MIN_THRESH
            continue;
        end
        
        [E,dE] = compute_E(dict{r}.A1,dict{id_atom}.A1);
        delta_L_C1 = delta_L_C1 + 2*tmpAlpha_rj*dict{id_atom}.C1*(E.*(dict{id_atom}.C1'*dict{r}.C1));
        delta_L_C1 = real(delta_L_C1);
        tmpC_rj = dict{r}.C1'*dict{id_atom}.C1;
        delta_L_A1 = delta_L_A1 + tmpAlpha_rj*sum(tmpC_rj.*dE'.*conj(tmpC_rj),2);
        delta_L_A1 = real(delta_L_A1 );
    end
    
    % compute gradients over data
    for id_data=1:nData
        if alpha(r,id_data)< MIN_THRESH
            continue;
        end
        [E_1, dE_1] = compute_E(dict{r}.A1, data{id_data}.A1);
%         [E_2, dE_2] = compute_E(dict{r}.A1, data{id_data}.A2);
        
%         dataQ_i2 = data{id_data}.C2 * (kron(eye(n/2),1/sqrt(2)*[1 -1i;1 1i])); % return to real domain to save compuation below
        tmpC_ri1 = dict{r}.C1'* data{id_data}.C1;
%         tmpC_ri2 = dict{r}.C1'* dataQ_i2; 
        
        delta_L_C1 = delta_L_C1 - 2*alpha(r,id_data)*data{id_data}.C1*(E_1.*(data{id_data}.C1'*dict{r}.C1));
%         delta_L_C1 = delta_L_C1 - 2*alpha(r,id_data)*dataQ_i2*(E_2.*(dataQ_i2'*dict{r}.C1));
        delta_L_A1 = delta_L_A1 - alpha(r,id_data)*sum(tmpC_ri1.*dE_1'.*conj(tmpC_ri1),2);
%         delta_L_A1 = delta_L_A1 - alpha(r,id_data)*sum(tmpC_ri2.*dE_2'.*conj(tmpC_ri2),2);
        
        delta_L_C1 = real(delta_L_C1);
        delta_L_A1 = real(delta_L_A1);
    end
    
    % normalize the gradients
    delta_L_C1 = 1/nData * delta_L_C1;
    delta_L_A1 = 1/nData * delta_L_A1;

    % update A1
    delta = SPGD_opts.learningRate_A;
    a = SPGD_opts.SN_a;
    b = SPGD_opts.SN_b;
    rho_old = 1/a*(log(b+dict{r}.A1)-log(b-dict{r}.A1));
    rho_new = rho_old -2*delta*a*b*delta_L_A1.*sigm(a*rho_old).*(1-sigm(a*rho_old));
    dict{r}.A1 = real(b*(2*sigm(a*rho_new)-1)); % in case that rho_old and rho_new are complex when previous lambda=1
    
    % update C1
    delta_t = SPGD_opts.learningRate_C;
    dict{r}.C1 = project_to_manifold(dict{r}.C1, delta_L_C1, delta_t);
end
  

newDict = dict;

end

