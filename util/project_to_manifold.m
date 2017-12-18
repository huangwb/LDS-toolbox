function Y = project_to_manifold(X, G, t)
%PROJECT_TO_MANIFOLD projects the gradients to the stefel manifold via Cayley transformation
%
% INUTS
% X               - [] initial point 
% G               - [] gradient in Euclidean space
% t               - [] learning rate
%
% OUTPUTS
% Y              - [] updated point for minimizing the objective
% implemented by Wenbing Huanng, 2016-6-12
p = size(G,2);

U = [G,X];
V = [X,-G];

Y = X - t*U/(eye(2*p)+t/2*V'*U)*V'*X;

end

