function d = distMartin( A1,C1,A2,C2)
%DISTMATIN: compute the Martin distance between two LDSs with
%system tuples (A1,C1) and (A2,C2)


O11 = dlyap(A1',C1'*C1);
O12 = dlyap(A1',A2,C1'*C2);
% O21 = dlyap(A2',A1,C2'*C1); %O21 = O12'
O22 = dlyap(A2',C2'*C2);

[~,S11,~] = svd(O11);
[~,S12,~] = svd(O12);
[~,S22,~] = svd(O22);

d = -sum(2*log(diag(S12)))+sum(log(diag(S11)))+sum(log(diag(S22)));
  

end

