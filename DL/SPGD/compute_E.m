function [E,dE] = compute_E(lambda_r,Lambda)
% lambda_r = r*a;
% E(j,i)  = (1-|lambda_r(i)|^2)(1-|Lambda(j)|^2)/|1-lambda_r(i)Lambda(j)^*|^2
% dE(j,.)/dr = [(1-|theta_j(1)|^2)/|1-theta_r(i)conj(theta_j(1))|^4*((1+r^2)(a*conj(Lambda(j))+conj(a)*Lambda(j))-2r(1+|Lambda|^2))

n = length(lambda_r);
lambda_r = repmat(lambda_r.',n,1);
Lambda = repmat(Lambda,1,n);

% r is real. lambda_r = r*a, if lambda_r is real, a = 1; if lambda_r is image, a=1i;
r = real(lambda_r)+imag(lambda_r);  
if isreal(lambda_r)
    a = 1;
else
    a = 1i;
end


tmp1 = 1-lambda_r.*conj(lambda_r);
tmp2 = 1-Lambda.*conj(Lambda);
tmp3 = (1-lambda_r.*conj(Lambda)).*(1-conj(lambda_r).*Lambda);
tmp4 = conj(a).*Lambda+a.*conj(Lambda);

E = (tmp1.*tmp2)./tmp3 ;
dE = tmp2./tmp3.^2.*((r.^2+1).*tmp4-2*r.*(2-tmp2));

E = real(E);
dE = real(dE);

end