function [LDS, scores, time] = dytex(Y,opts)
% Extracts dynamical testures from time series;
%
% INUTS
% Y           - {} input sequences stored in cell format 
% LDS_opts - .  LDS set-up parameters, see ./params/default.m
%
% OUTPUTS
% LDS        - .  LDS dynamical testures
% scores     - [] Returns the reconstruction errors (could be of multiple values recording the convergent process)
% time       - [] Returns the running time
%
% implemented by Wenbing Huanng, 2016-12-04

% MIN_THRESH = 0.01;

% Set-up
scores = [];
bound = opts.bound;
n = opts.dim_hidden;
nv = opts.nv;
svd_check = opts.svd_check;
tau = size(Y,2);

Y = Y + 1e-5*randn(size(Y)); % make sure Y is not equal to a zero matrix
Ymean = mean(Y,2);
Y = Y - Ymean*ones(1,tau); % Center input sequences
Y = Y';


% Perform SVD on Hankel Matrix
if size(Y,2) < size(Y,1)
    [V,S,U] = svd(Y,0);
else
    [U,S,V] = svd(Y',0);
end

% Estimatation of the State Vector X and the Measurement Matrix C
V = V(:,1:n);
S = S(1:n,1:n);
U = U(:,1:n);
Xhat = S*V';
Chat = U;

% Learn the transition matrix via different methods
begintime = clock; 
switch opts.transition
    
    case 'SYM' 
        [theta,P] = symTranMatrix(Xhat(:,1:tau-1),Xhat(:,2:tau));
        Ahat = theta;
        
        Vhat = (P*Xhat(:,2:tau)-diag(theta)*P*Xhat(:,1:tau-1));
        %     Vhat = diag(1./(diag(S(1:n,1:n))))*(P*Xhat(:,2:tau)-diag(theta)*P*Xhat(:,1:tau-1));
        [Uv,~,~] = svd(Vhat);
        Bhat = Uv(:,1:nv);
        Chat = Chat*P';
        
        
        switch opts.stabilizer
            case 'NONE'
                % no further computations     
            
            case 'SN'
                a = opts.SN_a;
                b = opts.SN_b;
                Ahat = b*(2*sigm(a*Ahat)-1);
                
            case 'UN'
                Ahat = 1/opts.UN_factor*Ahat;
        end

        LDS.A1 = Ahat;
%         LDS.B = Bhat;
        LDS.C1 = Chat;
%         LDS.R = Chat*Bhat;
        
    case 'FREE'
%         Ahat = S*V(2:end,:)'*V(1:end-1,:)*(eye(n)+1/(1-V(end,:)*V(end,:)')*V(end,:)'*V(end,:))*diag(1./diag(S));
        Ahat = Xhat(:,2:tau)*pinv(Xhat(:,1:tau-1));
        scores(end+1) = norm(S*V(2:end,:)'- Ahat*S*V(1:end-1,:)','fro')^2;
        
        switch opts.stabilizer
            case 'NONE'
                % no further computations     
            
            case 'SN'
                a = opts.SN_a;
                b = opts.SN_b;
                [U_A,S_A,V_A] = svd(Ahat);
                Ahat = U_A*diag(b*(2*sigm(a*diag(S_A))-1))*V_A';
                scores(end+1) = norm(S*V(2:end,:)'- Ahat*S*V(1:end-1,:)','fro')^2;
            
            case 'UN'
                Ahat = 1/opts.UN_factor*Ahat;
                
            case 'MN'
                Ahat = max(1/(1+(1-V(end,:)*V(end,:)')^(-0.5)),1-V(end,:)*V(end,:)')*Ahat;         
                scores(end+1) = norm(S*V(2:end,:)'- Ahat*S*V(1:end-1,:)','fro')^2;
                
            case 'HZP'
                Ahat = S*V(2:end,:)'*V(1:end-1,:)*diag(1./diag(S));      
                scores(end+1) = norm(S*V(2:end,:)'- Ahat*S*V(1:end-1,:)','fro')^2;

            case 'WLS'
                [Ahat, scores] = LearnWLSModel(V,S,bound,svd_check);
                
            case 'DWLS'
                [Ahat, scores] = LearnDWLSModel(V,S,bound,svd_check);
                
            case 'CG'
                [Ahat, scores] = learnCGModel(Xhat(:,1:end-1), Xhat(:,2:end),bound, svd_check);
                
            case 'LB2'
                [Ahat, scores] = learnLB2Model(Xhat(:,1:end-1), Xhat(:,2:end));
                
            case 'LB1'
                [Ahat, scores] = learnLB1Model(Xhat(:,1:end-1), Xhat(:,2:end));
        end
    
        % Compute B given A, C, Xhat
%         Vhat = (Xhat(:,2:tau)-Ahat*Xhat(:,1:tau-1));
%         [Uv,~,~] = svd(Vhat);
%         Bhat = Uv(:,1:nv);

        switch opts.decomposition
            
            case 'NONE'
                LDS.A = Ahat;
%                 LDS.B = Bhat;
                LDS.C = Chat;
                %         LDS.Q = cov(Vhat');
%                 LDS.R = Chat*Bhat;
                %         LDS.Ymean = Ymean;
                %         LDS.Xhat = Xhat;
                
                dot = dlyap(Ahat',eye(n));
                [Udot,Sdot,~] = svd(dot);
                invL = diag(1./sqrt(diag(Sdot)))*Udot';
                LDS.iL = invL;
                
            case 'CANONICAL'
                [U_A,S_A,V_A] = svd(Ahat);
                LDS.S = diag(S_A);
                LDS.V = V_A'*U_A;
                LDS.U = Chat*U_A;
                
                dot = dlyap((S_A*LDS.V)',eye(n));
                [Udot,Sdot,~] = svd(dot);
                invL = diag(1./sqrt(diag(Sdot)))*Udot';
                LDS.iL = invL;
            
            case 'SYMSKEW'
                A1 = 1/2*(Ahat + Ahat');
                A2 = 1/2*(Ahat - Ahat');
                
                [P1,D1] = eig(A1);
                [P2,D2] = eig(A2);
                
                C1 = Chat*P1;
                C2 = Chat*P2;
                
                LDS.A1 = diag(D1);
                LDS.A2 = diag(D2);
                LDS.C1 = C1;
                LDS.C2 = C2;
                
            case 'GRASS'
                LDS.O = Chat;
                for j=1:opts.Grass_k-1
                    LDS.O = [Chat; LDS.O*Ahat];
                end
                LDS.O = orth(LDS.O);
        end

        time = etime(clock,begintime);
        
end

end

