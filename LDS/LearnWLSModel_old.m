function [Ahat,scores] = LearnWLSModel(V,S,bound)
% LearnWLSModel leans the transition matrix of a Linear Dynamical System using
% the weighted-least-square method
%
% INPUTS
% V         - [tau, n] the right orthogonal matrix of the SVD for observation vectors
% S         - [n n] the scale matrix of the SVD for observation vectors
% bound     - [1] the upper bound of the spectral radius of the transition
%             matrix
%
% OUTPUTS
% Ahat      - [n n] the learned transition matrix
% scores    - records the objective values along the optimization track

% specify the experimental settings 
display = 0;
eps = 0;
dotprint = 1;
maxIter = 1000;
NumQP = 1000;


% we keep track of the following quantities over iterations,
% they can be plotted if desired ...
maxevals = [];
minevals = [];
scores = [];

begintime = clock; 

% first Ahat is learned unconstrained
fprintf('calculating initial A...\n');
[P,lambda,~] = svd(V(1:end-1,:)'*V(1:end-1,:));
Abar = P'*(V(2:end,:)'*V(1:end-1,:))*P*diag(1./diag(lambda));
Ahat = S*P*Abar*P'*diag(1./diag(S));
lsscore = norm(S*V(2:end,:)'- Ahat*S*V(1:end-1,:)','fro')^2;
scores(end+1) = lsscore;
fprintf('frob score = %.4f\n',lsscore);

% compute the eigenvalues of Abar instead of Ahat
[tmp_evals,max_e,min_e] = get_eigenthings(Abar);
maxevals(end+1) = max_e;
minevals(end+1) = min_e;

% do the initial check
if max_e < bound & min_e > -bound
    fprintf('done early!\n');
    fprintf('Top eigenval = %.4f\n',max_e);
    fprintf('Smallest eigenval = %.4f\n',min_e);
    return;
else
    fprintf('initial top eigenval = %.4f\n',max_e);
    fprintf('initial smallest eigenval = %.4f\n',min_e);
end


% calculate the required terms for QP: 
% W_H, W_f, W_A, W_b, W_c
constraints = 0;
M = S*(V(2:end,:)'*V(1:end-1,:))*P*diag(diag(lambda).^(-1/2));
iM = sum(M.^2);
W_H = diag(iM);
W_f = -iM';
W_A = [];
W_b = [];
W_c = trace(S*(V(2:end,:)'*V(2:end,:))*S); 
W = ones(size(V,2),1);
Abbar = Abar * diag(W);


% specify the QP optimization method with the interior-point-convex method
options = optimoptions('quadprog','Algorithm','interior-point-convex', 'Display','off','MaxIter',maxIter);
for i=1:NumQP
    Abbarprev = Abbar;
    [W_u,W_s,W_v] = svds(Abbar,1);
    tmp = diag(W_v*W_u'*Abar);
    W_A = [W_A; tmp'];
    W_b = [W_b; bound];
    [W,FVAL,EXITFLAG] = quadprog(W_H, W_f, W_A, W_b, [], [], [], [], [], options);
    Abbar = Abar * diag(W);
    
    % change the QP optimization method with the active-set method if the
    % interior-point-convex method fail to converge
    if EXITFLAG == -2
        instead_options = optimoptions('quadprog','Algorithm','active-set', 'Display','off','MaxIter',maxIter);
        [W,FVAL,EXITFLAG] = quadprog(W_H, W_f, W_A, W_b, [], [], [], [], [], instead_options);
        Abbar = Abar * diag(W);
    end
    
    score = (2*FVAL+W_c);
    scores(end+1) = score; 
    diffscore = (score - lsscore)/lsscore;
    
    [tmp_evals,max_e,min_e] = get_eigenthings(Abbar);
    
    maxevals = [maxevals max_e];
    minevals = [minevals min_e];
    
    if( max_e < bound & min_e > -bound )    % then we're done
        fprintf('found M, exiting ...\n');
        fprintf('eps: %.3f, diffscore: %.4f, top eval: %.7f, small eval: %.7f\n',eps,diffscore,max_e,min_e);
        fprintf('time so far: %.3f\n', etime(clock,begintime));
        break;
    else
        fprintf('.');
        if mod(i,dotprint) == 0 & display
            fprintf('\n');
            fprintf('eps: %.3f, diffscore: %.4f, top eval: %.7f, small eval: %.7f\n',eps,diffscore,max_e,min_e);
            fprintf('time so far: %.3f\n', etime(clock,begintime));
        end
    end
    constraints = i;
end

maxeig= max(abs(eig(Abbar)));
maxsval = svds(Abbar,1);

fprintf('top eigval after constrained QP: %.7f\n',maxeig);
fprintf('top sval after constrained QP: %.7f\n',maxsval);


% refining the solution: binary search to find boundary of stability region
Abbest = [];
tol = 0.00001;
lo = 0;
hi = 1;

Aborig = Abbarprev;
% Aborig = Abar;
while hi-lo > tol
    fprintf(',');
    alpha = lo + (hi-lo)/2;
    Abbest = (1-alpha)*Abbar + alpha*Aborig;
    maxeig = max(abs(eig(Abbest)));
    if (maxeig) > bound
        hi = alpha;
    elseif maxeig < bound
        lo = alpha;
    else    % done!
        break
    end
end
Abbest = (1-alpha + tol)*Abbar + (alpha-tol)*Aborig;
maxeig= max(abs(eig(Abbest)));
maxsval = svds(Abbest,1);

fprintf('After binary search: \n')

Ahat = S*P*Abbest*P'*diag(1./diag(S));
score = norm(S*V(2:end,:)'- Ahat*S*V(1:end-1,:)','fro')^2;
scores(end+1) = score;
diffscore = (score - lsscore)/lsscore;
fprintf('\n');
fprintf('eps: %.3f, diffscore: %.4f, top eval: %.7f, top sval: %.7f\n',eps,diffscore,maxeig,maxsval);
time = etime(clock,begintime);
fprintf('time so far: %.3f\n', time);

end

function [actual_evals,max_e,min_e] = get_eigenthings(A)
% maximum and minimum magnitudes of eigenvalues and corresponding
% eigenvectors, and the actual eigenvalues
[tmp_evecs,tmp_evals] = eig(A);
actual_evals = diag(tmp_evals);
evals = abs(actual_evals);
max_e = max(evals);
min_e = min(evals);

end
