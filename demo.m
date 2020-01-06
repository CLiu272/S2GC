clear
clc

addpath(genpath(pwd))

% please download the dataset first: https://pan.baidu.com/s/1JlnKE5ER9lsfqWBgL6fVog
load  BIC 
 
[Data1] = ProgressData(Gene);
[Data2] = ProgressData(Methy);
[Data3] = ProgressData(Mirna);

X = [Data1  Data2  Data3];
y = table2array(Response(:,2:end));
%y(:,1) = center(y(:,1));
X = [y X];

ViewNum(1) = size(Data1,2);
ViewNum(2) = size(Data2,2);
ViewNum(3) = size(Data3,2);

ViewIdx{1} = 1:ViewNum(1);
ViewIdx{2} = (ViewNum(1)+1) : (ViewNum(1)+ViewNum(2));
ViewIdx{3} = (ViewNum(1)+ViewNum(2)+1) : (ViewNum(1)+ViewNum(2)+ViewNum(3));

X = Cox_data_processed( X(:,3:end),X(:,1),'censoring',~X(:,2) );

lambda = 1;
rho = 0.1;
tau = 10;
[wm, funcValMv,G] = S2GC(X, lambda, rho, tau, ViewIdx);

OptimalClusterNum = 4;   
% BIC: 4 ; 
% COAD: 6 ;
% KRCCC:4 ;  7   
% LSCC: 3 ;   
% GBM: 3

idx = SpectralClustering(double(G),OptimalClusterNum);
group = num2str(idx);
group = num2cell(group);
[p,fh,stats]=MatSurv(X.sorty,~X.cens,group)



