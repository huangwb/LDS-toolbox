function cost = computeDictCost(data,dict,alpha,SC_opts)
%COMPUTEDICTCOST computes the objective values given data and dictionary

computeKernel = SC_opts.kernel;


K_D = computeKernel(dict);
K_XD = computeKernel(dict,data);


cost =  1/length(data)*(-2*trace(alpha'*K_XD) + trace(alpha'*K_D*alpha));

if cost>0
    fprintf('stop!');
%     save('debug.mat','data','dict','alpha');
%     return;
end

end