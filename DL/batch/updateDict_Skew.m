function newDict = updateDict_Skew(data,dict,alpha,LDS_opts)
%UPDATEDICT update the skew-symmetric part of the dictionary atoms
%
% INUTS
% data            - {} LDS data
% dict            - {} LDS dictionary
% alpha           - [] Sparse codes
%
% OUTPUTS
% newDict        - {} Updated dictionary
% implemented by Wenbing Huanng, 2016-4-23


MIN_THRESH = 1e-6;
nAtoms = length(dict);
[m,n] = size(dict{1}.C2);

%------ update A2 and C2 -------
for r=1:nAtoms
    
    alpha_nonzero = find(alpha(r,:),1);
    if (isempty(alpha_nonzero))
        fprintf('a useless atom identified!\n');
        continue;
    end
    for id_column=1:n/2
        
        Q_r = dict{r}.C2 * (kron(eye(n/2),1/sqrt(2)*[1,-1i;1,1i]));
           
        numFree = m-n;
        W = [null(Q_r'),Q_r(:,2*id_column-1),Q_r(:,2*id_column)] ;
        S = zeros(numFree+2);
        dS = 0;   
        for id_atom=1:nAtoms
            tmpAlpha_rj = alpha(r,:)*alpha(id_atom,:)';
            if id_atom == r || abs(tmpAlpha_rj) < MIN_THRESH
                continue;
            end
            
                       
            [E_odd, dE_odd] = compute_E_skew(dict{r}.A2(2*id_column-1), dict{id_atom}.A2);
            [E_even, dE_even] = compute_E_skew(dict{r}.A2(2*id_column), dict{id_atom}.A2);
%             E_old = (1-dict{r}.A2(id_column)^2)*((1-dict{id_atom}.A2.^2)./(1-dict{r}.A2(id_column)*dict{id_atom}.A2).^2) ;
%             dE_old = 2/(1-dict{r}.A2(id_column)^2)*(dict{id_atom}.A2-dict{r}.A2(id_column))./(1-dict{r}.A2(id_column)*dict{id_atom}.A2).*E;            
            
            
            tmpWC_j = W'*dict{id_atom}.C2;
            tmpC_rj_odd  = Q_r(:,2*id_column-1)'*dict{id_atom}.C2;
            tmpC_rj_even = Q_r(:,2*id_column)'*dict{id_atom}.C2;
            S = S + tmpAlpha_rj*1/2*(tmpWC_j.*(repmat(E_odd',[numFree+2,1]) + repmat(E_even',[numFree+2,1])))*tmpWC_j';
            dE = dE_odd - dE_even;
            dS = dS + tmpAlpha_rj*1/2*((tmpC_rj_odd .* dE')*tmpC_rj_odd'+ (tmpC_rj_even .* dE')*tmpC_rj_even');
        end
        for id_data=1:length(data)
            if alpha(r,id_data)< MIN_THRESH
                continue;
            end
 
%             [F_1_odd, dF_1_odd] = compute_E_skew(dict{r}.A2(2*id_column-1), data{id_data}.A1);
%             [F_1_even, dF_1_even] = compute_E_skew(dict{r}.A2(2*id_column), data{id_data}.A1);
            [F_2_odd, dF_2_odd] = compute_E_skew(dict{r}.A2(2*id_column-1), data{id_data}.A2);             
            [F_2_even, dF_2_even] = compute_E_skew(dict{r}.A2(2*id_column), data{id_data}.A2); 
                                       
%             F_old =  (1-dict{r}.A2(id_column)^2)*((1-data{id_data}.A2.^2)./(1-dict{r}.A2(id_column)*data{id_data}.A2).^2) ;
%             dF_old = 2/(1-dict{r}.A2(id_column)^2)*(data{id_data}.A2-dict{r}.A2(id_column))./(1-dict{r}.A2(id_column)*data{id_data}.A2).*F;                   
            
            dataQ_i2 = data{id_data}.C2 * (kron(eye(n/2),1/sqrt(2)*[1, -1i;1, 1i]));
%             tmpWC_i1 = W'*data{id_data}.C1;            
            tmpWC_i2 = W'*dataQ_i2;
%             tmpC_ri1_odd  = Q_r(:,2*id_column-1)'*data{id_data}.C1;
%             tmpC_ri1_even = Q_r(:,2*id_column)'*data{id_data}.C1;
            tmpC_ri2_odd  = Q_r(:,2*id_column-1)'*dataQ_i2;            
            tmpC_ri2_even = Q_r(:,2*id_column)'*dataQ_i2;           
            
%             S = S-alpha(r,id_data)*1/4*([tmpWC_i1, tmpWC_i2].*[(repmat(F_1_odd',[numFree+2,1]) + repmat(F_1_even',[numFree+2,1])),(repmat(F_2_odd',[numFree+2,1]) + repmat(F_2_even',[numFree+2,1]))]* [tmpWC_i1'; tmpWC_i2']);
%             dF_1 = dF_1_odd - dF_1_even;
%             dF_2 = dF_2_odd - dF_2_even;
%             dS = dS - alpha(r,id_data)*1/4*( (tmpC_ri1_odd.*dF_1')*tmpC_ri1_odd' + (tmpC_ri1_even.*dF_1')*tmpC_ri1_even'  + (tmpC_ri2_odd.*dF_2')*tmpC_ri2_odd' + (tmpC_ri2_even.*dF_2')*tmpC_ri2_even');
            
            S = S-alpha(r,id_data)*1/2*(tmpWC_i2.* (repmat(F_2_odd',[numFree+2,1]) + repmat(F_2_even',[numFree+2,1])) * tmpWC_i2');
            dF_2 = dF_2_odd - dF_2_even;
            dS = dS - alpha(r,id_data)*1/2*((tmpC_ri2_odd.*dF_2')*tmpC_ri2_odd' + (tmpC_ri2_even.*dF_2')*tmpC_ri2_even');            
        end
        % back to the real domain
        S = real(S);
        dS = real(dS);
%         u0 = zeros(numFree+1,1);
%         u0(end) = 1;
%         u = GCG( u0,S,20 );
        u = smallestEigvector(double(S),2);
        Wu = W*u;
        dict_Cr = Wu * 1/sqrt(2)*[1 1;1i -1i];
        dict{r}.C2(:,2*id_column-1) = dict_Cr(:,1);
        dict{r}.C2(:,2*id_column) = dict_Cr(:,2);

        lambda = real(dict{r}.A2(2*id_column-1)*(-1i));        
        delta = 0.0001;
        a = LDS_opts.SN_a;
        b = LDS_opts.SN_b;
        
        rho_old = 1/a*(log(b+lambda)-log(b-lambda));
        rho_new = rho_old -2*delta*dS*a*b*sigm(a*rho_old)*(1-sigm(a*rho_old));
        lambda = real(b*(2*sigm(a*rho_new)-1));       
        
        dict{r}.A2(2*id_column-1) = lambda*1i;
        dict{r}.A2(2*id_column)  = -lambda*1i;
        
        %    
    end
  
end
%----------------------------------------

newDict = dict;

end

function [E,dE] = compute_E_skew(theta_rk,theta_j)

lambda_rk = real(theta_rk*(-1i));  % lambda_rk is real. if A_r is symmetric, dict{r}.A2(k) = lambda_rk; if A_r is skew-symmetric, dict{r}.A2(k) = lambda_rk*i

tmp1 = 1-theta_rk.*conj(theta_rk);
tmp2 = 1-theta_j.*conj(theta_j);
tmp3 = (1-theta_rk*conj(theta_j)).*(1-conj(theta_rk)*theta_j);
tmp4 = theta_j*1i+conj(theta_j*1i);

E = tmp1*(tmp2./tmp3) ;
dE = tmp2./tmp3.^2.*((lambda_rk^2+1)*tmp4-2*lambda_rk*(2-tmp2));

E = real(E);
dE = real(dE);
end
