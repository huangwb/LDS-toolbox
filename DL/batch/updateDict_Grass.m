function newDict = updateDict_Grass(data,dict,alpha,LDS_opts)



MIN_THRESH = 1e-6;
nAtoms = length(dict);
[mt,n] = size(dict{1}.O);

%------ update A1 and C1 -------
for r=1:nAtoms
    S = zeros(mt);
    idx_alpha = find(alpha(r,:));
    if (isempty(idx_alpha))
        fprintf('a useless atom identified!\n');
        continue;
    end
    
    for tmpC1 = 1:length(idx_alpha)
        S = S -alpha(r,idx_alpha(tmpC1))*data{idx_alpha(tmpC1)}.O*data{idx_alpha(tmpC1)}.O';
        for j = 1:nAtoms
            if (j == r) || (alpha(j,idx_alpha(tmpC1)) == 0)
                continue;
            end
            S = S + alpha(j,idx_alpha(tmpC1))*alpha(r,idx_alpha(tmpC1))*dict{j}.O*dict{j}.O';
        end
    end
    
    [O,~] = eigs(double(S),n);
    dict{r}.O = O;
    
end


newDict = dict;

end



