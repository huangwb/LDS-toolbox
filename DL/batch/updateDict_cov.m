function newDict = updateDict_cov(data,dict,alpha,LDS_options)
%updateDict_cov updates the covariances of the dictioanry
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
[m,n] = size(dict{1}.C);
nv = size(dict{1}.R,2);


%------ update R (project to B) -------
for r=1:nAtoms
    S = zeros(m);
    for id_atom=1:nAtoms
        tmpAlpha_rj = alpha(r,:)*alpha(id_atom,:)';
        if id_atom == r || abs(tmpAlpha_rj) < MIN_THRESH
            continue;
        end
        S =  S + tmpAlpha_rj * dict{id_atom}.R * dict{id_atom}.R';  %#ok<MHERM>
        
    end
    
    for id_data=1:length(data)
        if alpha(r,id_data)==0
            continue;
        end
        S = S-alpha(r,id_data)* data{id_data}.R * data{id_data}.R';   %#ok<MHERM>       
    end
    
   R = smallestEigvector(S,nv);
   dict{r}.R = R;
    
end

newDict = dict;

end

