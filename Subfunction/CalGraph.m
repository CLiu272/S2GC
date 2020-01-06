function [L,SS] = CalGraph(Data,w,k,tau)

mu = 10e-1;
alpha = 1;
n = size(Data{1}.X,1);
ViewNum = length(Data);
distX = zeros(n,n);
distf = zeros(n,n);

for m = 1:ViewNum
    distX =  distX + tau(m) * L2_distance_1(Data{m}.X',Data{m}.X');    
    y = Data{m}.X*w{m};
    distf = distf + tau(m) * L2_distance_1(y',y');
end


[~, idx] = sort(mu*distX+distf,2);

S = zeros(n,n);

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

