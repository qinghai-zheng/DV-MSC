function[S] = msc(X,opts)
addpath('./LRR')
addpath('./twist')
%% preparation...

V = size(X,2);
N = size(X{1},2);

lambda_1 =opts.lambda_1;
lambda_2 =opts.lambda_2;
d_k = opts.dim_k;

P = cell(1,V);
E1 = cell(1,V);
E2 = cell(1,V);
Q = cell(1,V);
Z = cell(1,V);
Y1 = cell(1,V);
Y2 = cell(1,V);
Y3 = cell(1,V);
H = cell(1,V);
K = cell(1,V);
dim_all_view = zeros(V,1);

SD = 0;

for i=1:V
    dim_all_view(i) = size(X{i},1);
    SD = SD + dim_all_view(i);
    P{i} = zeros(dim_all_view(i),d_k);
    E1{i} = zeros(dim_all_view(i),N);
    E2{i} = zeros(d_k,N);
    Q{i} = zeros(N,N);
    Z{i} = zeros(N,N);
    Y1{i} = zeros(dim_all_view(i),N);
    Y2{i} = zeros(d_k,N);
    Y3{i} = zeros(N,N);
    H{i} = zeros(d_k,N);
end

K = updateK(H,V,N);

mu = 10e-6; max_mu = 10e6;

eta = 1.1;

isConverge = 0;
epsilon = 1e-6;

sX = [N, N, V];
iter = 0;
max_iter = 250;

%% updating variables...
while (isConverge == 0&&iter<=max_iter)

    for i=1:V
        % update P{i}
        
        P{i} = updateP(Y1{i},mu,X{i},H{i},E1{i});
        
        % update H{i}
        
        H{i} = updateH(X{i},P{i},Z{i},E1{i},E2{i},Y1{i},Y2{i},K{i},mu,mu,lambda_1,N);
        
        % update Z{i}
        
        Z{i} = (mu*H{i}'*H{i}+mu*eye(N))\(H{i}'*Y2{i}-mu*H{i}'*E2{i}+mu*H{i}'*H{i}-Y3{i}+mu*Q{i});
        
        % update E{i}
        a = X{i}-P{i}*H{i}+Y1{i}/mu;
        b = H{i}-H{i}*Z{i}+Y2{i}/mu;
        T_E = [X{i}-P{i}*H{i}+Y1{i}/mu; H{i}-H{i}*Z{i}+Y2{i}/mu];
        [E] = solve_l1l2(T_E,1/mu);
        E1{i} = E(1:dim_all_view(i),:);
        E2{i} = E(1+dim_all_view(i):dim_all_view(i)+d_k,:);
        
        % update Y1{i}
        
        Y1{i} = Y1{i} + mu*(X{i}-P{i}*H{i}-E1{i});
        
        % update Y2{i}
        
        Y2{i} = Y2{i} + mu*(H{i}-H{i}*Z{i}-E2{i});
        
    end
    
    % update K
    
    K = updateK(H,V,N);
    
    % update Q
    
    Z_tensor = cat(3, Z{:,:});
    Y3_tensor = cat(3, Y3{:,:});
    z = Z_tensor(:);
    y3 = Y3_tensor(:);
    
    %twist-version
    [q, objV] = wshrinkObj(z + 1/mu*y3,lambda_2/mu,sX,0,3);
    Q_tensor = reshape(q, sX);
    
    %update y3
    
    y3 = y3 + mu*(z - q);
    
    %record the iteration information
    history.objval(iter+1)   =  objV;
    
    %% coverge condition
    isConverge = 1;
    for k=1:V
        if (norm(X{k}-P{k}*H{k}-E1{k},inf)>epsilon)
            history.norm_H = norm(X{k}-P{k}*H{k}-E1{k},inf);
%                         fprintf('    norm_H %7.10f    ', history.norm_H);
            isConverge = 0;
        end
        
        if (norm(H{k}-H{k}*Z{k}-E2{k},inf)>epsilon)
            history.norm_Z = norm(H{k}-H{k}*Z{k}-E2{k},inf);
%                         fprintf('    norm_Z %7.10f    ', history.norm_Z);
            isConverge = 0;
        end
        
        Q{k} = Q_tensor(:,:,k);
        Y3_tensor = reshape(y3, sX);
        Y3{k} = Y3_tensor(:,:,k);
        if (norm(Z{k}-Q{k},inf)>epsilon)
            history.norm_Z_Q = norm(Z{k}-Q{k},inf);
%                         fprintf('norm_Z_G %7.10f    \n', history.norm_Z_Q);
            isConverge = 0;
        end
    end
    
    iter = iter + 1;
    mu = min(mu*eta, max_mu);
    
end
S = 0;
for k=1:V
    S = S + abs(Z{k})+abs(Z{k}');
end


