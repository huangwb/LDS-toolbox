function [M, scores] = learnLB2Model(S1,S2)
% Implements algorithm in
% Seth L. Lacy and Dennis S. Bernstein. Subspace Identification with
% Guaranteed Stability using Constrained Optimization. IEEE Transactions on
% Automatic Control, 48(7):1259-1263, July 2003
%
% 
% CVX (with Sedumi), publicly available
%   optimization software, is used to solve the optimization problem
%
% Authors: Sajid Siddiqi and Byron Boots
begintime = clock; 

scores = [];
n = size(S1,1);
L = size(S1,2);
r1 = n;
s = n;
delta = .001;
lambda = 0;


% we keep track of the following quantities over iterations,
% they can be plotted if desired ...
maxevals = [];
minevals = [];
scores = [];
svals = [];


% first M is learned unconstrained
fprintf('calculating initial M...\n');
M = pinv(S1')*S2';
lsscore = norm(S1'*M - S2','fro')^2;
fprintf('frob score = %.4f\n',lsscore);


[tmp_evals,max_e,min_e] = get_eigenthings(M);
maxevals(end+1) = max_e;
minevals(end+1) = min_e;
scores(end+1) = lsscore;

% watch for stability
if max_e < 1 & min_e > -1
    fprintf('done early! full eigenvals:\n');
    M = M';
    return;
else
    fprintf('initial top eigenval = %.4f\n',max_e);
    fprintf('initial smallest eigenval = %.4f\n',min_e);
end



U = zeros(size(S1));

tmp1 = [S1 ; U];
term2 = tmp1'*pinv(tmp1*tmp1');
term1 = S2;
term3 = term1*term2;

X1 = term3(1:n,1:n);
B = term3(1:n,n+1:end);

t1 = zeros(r1*s,1);
t2 = -eye(r1*s);
t3 = kron( [zeros(r1,n) eye(n)],[-eye(n) eye(n)*X1] );

t4 = zeros(n*n,1);
t5 = zeros(n*n,r1*s);
t6 = kron([zeros(n) eye(n)],[zeros(n) eye(n)]) - kron([eye(n) zeros(n)],[eye(n) zeros(n)]);
Ax = [t1 t2 t3   ; ...
      t4 t5 t6];
  
bx = delta*[zeros(r1*s,1) ; reshape(eye(n), numel(eye(n)), 1)];

N = 4*n*n + r1*s + 1;

z1_inds = 1;
z2_inds = (1+1):(s*r1 + 1);
z3_inds = (s*r1+1 + 1): (s*r1+1 + 4*n*n);
tmp = reshape((s*r1+1 + 1):(s*r1+1 + 4*n*n),2*n,2*n);
z3_P_inds = tmp(n+1:2*n,n+1:2*n);
z3_P_diag_inds = diag(tmp(n+1:2*n,n+1:2*n));

cx = zeros(N,1);
cx(1) = 1;
cx(z3_P_diag_inds) = lambda;

cvx_begin

  variable z(N)
  minimize (cx'*z)
  subject to
     
     Ax*z == bx;
     z(z1_inds) >= norm(z(z2_inds));
     reshape(z(z3_inds),2*n,2*n) == semidefinite(2*n);
     
cvx_end

Z3 = reshape(z(z3_inds),2*n,2*n);

P = Z3(n+1:2*n,n+1:2*n);
Q = Z3(1:n,n+1:2*n);
Ahat = Q*inv(P);

M = Ahat;

maxeig= max(abs(eig(M)));
maxsval = svds(M,1);

fprintf('top eigval after CVX: %.7f\n',maxeig);
fprintf('top sval after CVX: %.7f\n',maxsval);

score = norm(S2- M*S1,'fro')^2;
scores(end+1) = score;
diffscore = (score - lsscore)/lsscore;
fprintf('\n');
fprintf('eps: %.3f, diffscore: %.4f, top eval: %.7f, top sval: %.7f\n',eps,diffscore,maxeig,maxsval);
fprintf('time so far: %.3f\n', etime(clock,begintime));


end

function [actual_evals,max_e,min_e] = get_eigenthings(M)
% maximum and minimum magnitudes of eigenvalues and corresponding
% eigenvectors, and the actual eigenvalues
[tmp_evecs,tmp_evals] = eig(M);
actual_evals = diag(tmp_evals);
evals = abs(actual_evals);
max_e = max(evals);
min_e = min(evals);

end
