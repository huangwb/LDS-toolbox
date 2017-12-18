function newDict = update_Skew(data,dict,alpha,SPGD_opts)
%SGD_SYM Update the skew-symmetric dictionary via SGD
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
n = size(dict{1}.A2,1);
%------ update A2 and C2 -------
for r=1:nAtoms
    
    alpha_nonzero = find(alpha(r,:),1);
    if (isempty(alpha_nonzero))
        fprintf('a useless atom identified!\n');
        continue;
    end
    delta_L_Q = 0;
    delta_L_A2 = 0;
    Q = dict{r}.C2 * (kron(eye(n/2),1/sqrt(2)*[1,-1i;1,1i]));
    swap_Q = zeros(size(Q));
    % swap the odd colunms with the even ones
    swap_Q(:,1:2:n-1) = Q(:,2:2:n);
    swap_Q(:,2:2:n) = Q(:,1:2:n-1);
    
    % compute gradients over dictionary
    for id_atom=1:nAtoms
        tmpAlpha_rj = alpha(r,:)*alpha(id_atom,:)';
        if id_atom == r || abs(tmpAlpha_rj) < MIN_THRESH
            continue;
        end
        
        [E,dE] = compute_E(dict{r}.A2,dict{id_atom}.A2);
        
        % compute gradients of S_r w.r.t. Q
        E_S = kron(E(:,1:2:n-1)+E(:,2:2:n),[1,1]);
        delta_L_Q = delta_L_Q + tmpAlpha_rj*dict{id_atom}.C2*(E_S.*(dict{id_atom}.C2'*Q));           
        % compute gradients of sigma_r w.r.t. Q
        E_sigma = kron(E(:,1:2:n-1)-E(:,2:2:n),[1,-1]);
        delta_L_Q = delta_L_Q + tmpAlpha_rj*dict{id_atom}.C2*(E_sigma.*(dict{id_atom}.C2'*swap_Q))*1i;

        
        % compute gradients w.r.t. A2
        dE(:,2:2:n)=-dE(:,2:2:n); % note that the odd are conjugate to the even
        tmpC_rj = dict{r}.C2'*dict{id_atom}.C2;
        delta_L_A2 = delta_L_A2 + tmpAlpha_rj*sum(tmpC_rj.*dE'.*conj(tmpC_rj),2);
        delta_L_A2 = real(delta_L_A2);
    end
    
    % compute gradients over data
    for id_data=1:nData
        if alpha(r,id_data)< MIN_THRESH
            continue;
        end
%         [E_1, dE_1] = compute_E(dict{r}.A2, data{id_data}.A1);
        [E_2, dE_2] = compute_E(dict{r}.A2, data{id_data}.A2); 
        dataQ_i2 = data{id_data}.C2 * (kron(eye(n/2),1/sqrt(2)*[1 -1i;1 1i])); % return to real domain to save compuation below
        
        % compute gradients of S_r w.r.t. Q
%         E_1_S = kron(E_1(:,1:2:n-1)+E_1(:,2:2:n),[1,1]);
%         delta_L_Q = delta_L_Q - alpha(r,id_data)*data{id_data}.C1*(E_1_S.*(data{id_data}.C1'*Q));
        E_2_S = kron(E_2(:,1:2:n-1)+E_2(:,2:2:n),[1,1]);
        delta_L_Q = delta_L_Q - alpha(r,id_data)*dataQ_i2*(E_2_S.*(dataQ_i2'*Q));
        % compute gradients of sigma_r w.r.t. Q
        E_2_sigma = kron(E_2(:,1:2:n-1)-E_2(:,2:2:n),[1,-1]);
        delta_L_Q = delta_L_Q + tmpAlpha_rj*data{id_data}.C2*(E_2_sigma.*(data{id_data}.C2'*swap_Q))*1i;
        
        % compute gradients w.r.t. A2
%         dE_1(:,2:2:n)=-dE_1(:,2:2:n); % note that the odd are conjugate to the even
        dE_2(:,2:2:n)=-dE_2(:,2:2:n); % note that the odd are conjugate to the even
%         tmpC_ri1 = dict{r}.C2'* data{id_data}.C1;
        tmpC_ri2 = dict{r}.C2'* data{id_data}.C2; 
%         delta_L_A2 = delta_L_A2 - alpha(r,id_data)*sum(tmpC_ri1.*dE_1'.*conj(tmpC_ri1),2);
        delta_L_A2 = delta_L_A2 - alpha(r,id_data)*sum(tmpC_ri2.*dE_2'.*conj(tmpC_ri2),2);
        delta_L_A2 = real(delta_L_A2);
    end
    
    % normalize the gradients
    delta_L_Q = 1/nData * delta_L_Q;
    delta_L_A2 = 1/nData * delta_L_A2;

    % update A2
    lambda = real(dict{r}.A2(1:2:n-1)*(-1i));
    delta = SPGD_opts.learningRate_A; 
    a = SPGD_opts.SN_a;
    b = SPGD_opts.SN_b;
    rho_old = 1/a*(log(b+lambda)-log(b-lambda));
    rho_new = rho_old -2*delta*a*b*((delta_L_A2(1:2:n-1)+delta_L_A2(2:2:n))).*sigm(a*rho_old).*(1-sigm(a*rho_old));
    lambda = real(b*(2*sigm(a*rho_new)-1)); % in case that rho_old and rho_new are complex when previous lambda=1
    dict{r}.A2(1:2:n-1) = lambda*1i;
    dict{r}.A2(2:2:n) = -lambda*1i;
    
    % update C2
    delta_t = SPGD_opts.learningRate_C;
    dict{r}.C2 = project_to_manifold(Q, delta_L_Q, delta_t)*kron(eye(n/2),1/sqrt(2)*[1,1;1i,-1i]);
end
  

newDict = dict;

end


