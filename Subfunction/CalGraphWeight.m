function [L,SS] = CalGraphWeight(X,w,k)

mu = 10e-4;
alpha = 1;
n = size(X,1);

distX =  L2_distance_1(X',X');

y = X*w;
distf = L2_distance_1(y',y');
S = zeros(n);

[~, idx] = sort(mu*distX+distf,2);

for i=1:n
    idxa0 = idx(i,2:k+1);
    dfi = distf(i,idxa0);
    dxi = distX(i,idxa0);
    distK = sum((mu*dxi+dfi)/(alpha))/k;
    distFull = mu*distX(i,:) + distf(i,:);
    d = (distK - distFull);
    d(find(d<=0)) = 0;
    d(i) = 0;
    S(i,:) = d;
    S(i,:) = d./sum(d);    
end

S = single(S);
SS = (S+S')/2;
D = diag(sum(SS));
L = D-SS;
        
end

