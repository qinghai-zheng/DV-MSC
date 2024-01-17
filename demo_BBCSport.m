clear;clc;
addpath('./ClusteringMeasure');
load('./dataset/my_BBCSport.mat');
fprintf('msc on BBCSport dataset\n');
numC = size(unique(gt),1);

X{1} = X1;
X{2} = X2;
for i = 1:2 
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
end

NMI_all = [];
ACC_all = [];
F_all = [];
AVG_all = [];
P_all = [];
RI_all = [];
opts.lambda_1 = 1; 
opts.lambda_2 = 1;
opts.dim_k = 100;

fprintf('lambda_1 = %f,lambda_2 = %f,dim_k = %f\n', opts.lambda_1,opts.lambda_2,opts.dim_k);

for i = 1:30
    fprintf('Clustering in %d-th iteration\n',i);
    S  = msc(X,opts);
    [NMI,ACC,F,AVG,P,RI]=clustering(abs(S)+abs(S'), numC, gt);
    fprintf('\tNMI: %f, ACC: %f, F: %f, AVG: %f, P: %f, RI: %f\n',NMI,ACC,F,AVG,P,RI);
    NMI_all = [NMI_all, NMI];
    ACC_all = [ACC_all, ACC];
    F_all = [F_all, F];
    AVG_all = [AVG_all, AVG];
    P_all = [P_all, P];
    RI_all = [RI_all, RI];
end
fprintf('---------------Average Results--------------\n');
fprintf('--------------------------------------------\n');
fprintf('NMI: %f(%f), ACC: %f(%f), F: %f(%f), AVG: %f(%f), P: %f(%f), RI: %f(%f)\n',mean(NMI_all),std(NMI_all),mean(ACC_all),std(ACC_all),mean(F_all),std(F_all),mean(AVG_all),std(AVG_all),mean(P_all),std(P_all),mean(RI_all),std(RI_all));
fprintf('--------------------------------------------\n');
fprintf('--------------------------------------------\n');
 
