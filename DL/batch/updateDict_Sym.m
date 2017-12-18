function newDict = updateDict_Sym(data,dict,alpha,LDS_opts)
%UPDATEDICT update the symmetric part of the dictionary atoms
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
[m,n] = size(dict{1}.C1);

%------ update A1 and C1 -------
for r=1:nAtoms
    
    alpha_nonzero = find(alpha(r,:),1);
    if (isempty(alpha_nonzero))
        fprintf('a useless atom identified!\n');
        continue;
    end
    for id_column=1:n
        
        numFree = m-n;
        W = [null(dict{r}.C1') dict{r}.C1(:,id_column)] ;
        S = zeros(numFree+1);
        dS = 0;
                
        for id_atom=1:nAtoms
            tmpAlpha_rj = alpha(r,:)*alpha(id_atom,:)';
            if id_atom == r || abs(tmpAlpha_rj) < MIN_THRESH
                continue;
            end
            
            [E, dE] = compute_E_sym(dict{r}.A1(id_column),dict{id_atom}.A1);
%             E_old = (1-dict{r}.A1(id_column)^2)*((1-dict{id_atom}.A1.^2)./(1-dict{r}.A1(id_column)*dict{id_atom}.A1).^2) ;
%             dE_old = 2/(1-dict{r}.A1(id_column)^2)*(dict{id_atom}.A1-dict{r}.A1(id_column))./(1-dict{r}.A1(id_column)*dict{id_atom}.A1).*E;            
            
            tmpWC_j = W'*dict{id_atom}.C1;
            tmpC_rj = dict{r}.C1(:,id_column)'*dict{id_atom}.C1;
            S = S + tmpAlpha_rj*1/2*(tmpWC_j.*repmat(E',[numFree+1,1]))*tmpWC_j';
            dS = dS + tmpAlpha_rj*1/2*(tmpC_rj .*dE')*tmpC_rj';
 
        end
        for id_data=1:length(data)
            if alpha(r,id_data)< MIN_THRESH
                continue;
            end
            
            [F_1, dF_1] = compute_E_sym(dict{r}.A1(id_column), data{id_data}.A1);
%             [F_2, dF_2] = compute_E_sym(dict{r}.A1(id_column), data{id_data}.A2);                            
%             F1_old =  (1-dict{r}.A1(id_column)^2)*((1-data{id_data}.A1.^2)./(1-dict{r}.A1(id_column)*data{id_data}.A1).^2) ;
%             dF1_old = 2/(1-dict{r}.A1(id_column)^2)*(data{id_data}.A1-dict{r}.A1(id_column))./(1-dict{r}.A1(id_column)*data{id_data}.A1).*F1_old;                   
            
%             dataQ_i2 = data{id_data}.C2 * (kron(eye(n/2),1/sqrt(2)*[1 -1i;1 1i]));
            tmpWC_i1 = W'*data{id_data}.C1;
            tmpC_ri1 = dict{r}.C1(:,id_column)'*data{id_data}.C1;
%             tmpWC_i2 = W'*dataQ_i2;
%             tmpC_ri2 = dict{r}.C1(:,id_column)'* dataQ_i2; 
            
%             test1 = alpha(r,id_data)*1/2*( [tmpWC_i1 tmpWC_i2].*[repmat(F_1',[numFree+1,1]), repmat(F_2',[numFree+1,1])] )* [tmpWC_i1'; tmpWC_i2'];
%             tmpWC_i2 = W'* data{id_data}.C2;
%             test2 = alpha(r,id_data)*1/2*( [tmpWC_i1 tmpWC_i2].*[repmat(F_1',[numFree+1,1]), repmat(F_2',[numFree+1,1])] )* [tmpWC_i1'; tmpWC_i2'];
%             S = S-alpha(r,id_data)*1/2*( [tmpWC_i1, tmpWC_i2].*[repmat(F_1',[numFree+1,1]), repmat(F_2',[numFree+1,1])] )* [tmpWC_i1'; tmpWC_i2'];
%             dS = dS-alpha(r,id_data)*1/2*([tmpC_ri1,tmpC_ri2].*[dF_1', dF_2'])* [tmpC_ri1'; tmpC_ri2'];
            S = S-alpha(r,id_data)*1/2*( tmpWC_i1.* repmat(F_1',[numFree+1,1]) )* tmpWC_i1';
            dS = dS-alpha(r,id_data)*1/2*( tmpC_ri1.* dF_1')* tmpC_ri1';            
            
        end
        % back to the real domain
        S = real(S);
        dS = real(dS);
%         u0 = zeros(numFree+1,1);
%         u0(end) = 1;
%         u = GCG( u0,S,20 );
        u = smallestEigvector(double(S),1);
        dict{r}.C1(:,id_column) = W*u;
            
        delta = 0.0001;
        a = LDS_opts.SN_a;
        b = LDS_opts.SN_b;
        rho_old = 1/a*(log(b+dict{r}.A1(id_column))-log(b-dict{r}.A1(id_column)));
        rho_new = rho_old -2*delta*dS*a*b*sigm(a*rho_old)*(1-sigm(a*rho_old));
        dict{r}.A1(id_column) = real(b*(2*sigm(a*rho_new)-1));
        %    
    end
  
end


newDict = dict;

end

function [E,dE] = compute_E_sym(theta_rk,theta_j)

lambda_rk = theta_rk;  % lambda_rk is real. if A_r is symmetric, dict{r}.A1(k) = lambda_rk; if A_r is skew-symmetric, dict{r}.A1(k) = lambda_rk*i

tmp1 = 1-theta_rk.*conj(theta_rk);
tmp2 = 1-theta_j.*conj(theta_j);
tmp3 = (1-theta_rk*conj(theta_j)).*(1-conj(theta_rk)*theta_j);
tmp4 = theta_j+conj(theta_j);

E = tmp1*(tmp2./tmp3) ;
dE = tmp2./tmp3.^2.*((lambda_rk^2+1)*tmp4-2*lambda_rk*(2-tmp2));

E = real(E);
dE = real(dE);
end

