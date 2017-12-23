% refresh the paramenters related to the newly-set ones

switch LDS_opts.decomposition
    case 'SYMSKEW'
        LDS_opts.kernel = @computeKernel_SymSkew;
    case 'CANONICAL'
        LDS_opts.kernel = @computeKernel_CanProject;
    case 'GRASS'
        LDS_opts.kernel = @computeKernel_Grass;
end
if LDS_opts.beta ~= 1
    LDS_opts.kernel = combine_kernel(SC_opts.kernel, @computeKernel_Cov);
    % Compute the kernel between state covariances. See definition in ./LDS/ 
end

SC_opts.kernel = LDS_opts.kernel;

Clustering_opts.kernel = LDS_opts.kernel;

Classifier_opts.SC_opts = SC_opts;

