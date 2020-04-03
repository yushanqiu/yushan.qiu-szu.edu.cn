function [disdis,entropy,Distance_norm,s]=dimselect(types)
switch types
    case 'g'
load gene_EMT;
A = gene_EMT';
    case 'r'
load RBP_EMT;
A = data';
end


nL = zeros(304,1);
nL(1:160)=1;
N = 304;
Dis = squareform(pdist(A));
H = eye(N)-ones(N)/N;
K = H*A*A'*H;
[Distance_norm,entropy,s]=entropyss(K);
inds=find(entropy==min(entropy));
disdis = Distance_norm(1:end-1)-Distance_norm(2:end);

