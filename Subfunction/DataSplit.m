function [ Xtrain,Xtest ] = DataSplit( X, percent )
%DATASPLIT 

if percent > 1 || percent <0
    error('splitting percentage error')
end


task_sample_size = size(X,1);
tSelIdx = randperm(task_sample_size) < task_sample_size * percent;
    
Xtrain = X(tSelIdx,:);
Xtest = X(~tSelIdx,:);

Xtrain = Cox_data_processed( Xtrain(:,3:end),Xtrain(:,1),'censoring',~Xtrain(:,2) );
    

end

