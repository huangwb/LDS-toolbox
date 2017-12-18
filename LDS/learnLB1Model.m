function [M scores] = learnLB1Model(S1,S2)
% Implements algorithm in
% Seth L. Lacy and Dennis S. Bernstein. Subspace Identification with
% Guaranteed Stability using Constrained Optimization. In Proc American
% Control Conference, 2002
%
% CVX (with Sedumi), publicly available
%   optimization software, is used to solve the optimization problem
%
% Authors: Sajid Siddiqi and Byron Boots
begintime = clock; 

n = size(S1,1);
L = size(S1,2);
r1 = n;
s = n;
m = n;

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

SU = [S1 ; U];
SU_R = SU'*pinv(SU*SU');

t1 = zeros(n*L,1);
t2 = kron(eye(n),piperp(SU))*permtrans(L,n);
t3 = zeros(n*L,n*n);
t4 = t3;
t5 = zeros(n*L,4*n*n);
t6 = zeros(n*n,1);
t7 = permtrans(n,n)*kron(eye(n),(SU_R*[eye(n) ; zeros(m,n)])')*permtrans(L,n);
t8 = zeros(n*n);
t9 = t8;
t10 = [zeros(n*n) eye(n*n)]*kron(eye(2*n),[eye(n) zeros(n)]);
t11 = t6;
t12 = t7;
t13 = t8;
t14 = t8;
t15 = permtrans(n,n)*[eye(n*n) zeros(n*n)]*kron(eye(2*n),[zeros(n) eye(n)]);
t16 = t6;
t17 = zeros(n*n,L*n);
t18 = eye(n*n);
t19 = t8;
t20 = -[eye(n*n) zeros(n*n)]*kron(eye(2*n),[eye(n) zeros(n)]);
t21 = t6;
t22 = t17;
t23 = t18;
t24 = t18;
t25 = zeros(n*n,4*n*n);
t26 = t6;
t27 = t17;
t28 = t8;
t29 = t8;
t30 = [zeros(n*n) eye(n*n)]*kron(eye(2*n),[zeros(n) eye(n)]);

Ax = [t1 t2 t3 t4 t5;
     t6 t7 t8 t9 t10;
     t11 t12 t13 t14 t15;
     t16 t17 t18 t19 t20;
     t21 t22 t23 t24 t25;
     t26 t27 t28 t29 t30];
 
% tt1 = vec(piperp(SU)*S2');
% tt2 = vec(S2*SU_R*[ eye(n) ; zeros(m,n) ])
tmp1 = piperp(SU)*S2';
tt1 = tmp1(:);
tmp2 = S2*SU_R*[ eye(n) ; zeros(m,n) ];
tt2 = tmp2(:);

tt3 = tt2;
tt4 = zeros(n*n,1);
tmp5 = eye(n); 
tt5 = tmp5(:);
tt6 = tt5;

bx = [ tt1 ; tt2 ; tt3 ; tt4 ; tt5 ; tt6 ];

N = 6*n*n + n*L + 1;
cx = zeros(N,1);
cx(1) = 1;

z1_inds = 1;
z2_inds = (z1_inds(end) + 1):(z1_inds(end) + n*L);
z3_inds = (z2_inds(end) + 1): (z2_inds(end) + n*n);
z4_inds = (z3_inds(end) + 1): (z3_inds(end) + n*n);
z5_inds = (z4_inds(end) + 1): (z4_inds(end) + 4*n*n);

cvx_begin

  variable z(N)
  minimize (cx'*z)
  subject to
     
     Ax*z == bx;
     z(z1_inds) >= norm(z(z2_inds));
     reshape(z(z3_inds),n,n) == semidefinite(n);    
     reshape(z(z4_inds),n,n) == semidefinite(n);    
     reshape(z(z5_inds),2*n,2*n) == semidefinite(2*n);    
cvx_end

Z5 = reshape(z(z5_inds),2*n,2*n);

Ahat = Z5(1:n,n+1:2*n);

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



function X = piperp(M)
    n = size(M,2);
    X = eye(n) - M'*pinv(M*M')*M;
end    



function Pnn = permtrans(m,n)
    Pnn = zeros(m*n,m*n);
    inds = reshape(1:m*n,n,m)';
%     inds = vec(inds);
    inds = inds(:);
    for i = 1:m*n
        Pnn(i,inds(i)) = 1;
    end
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
