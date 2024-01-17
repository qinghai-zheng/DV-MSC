function P=updateP(Y1_i,mu,X_i,H_i,E1_i)
A = H_i';
B = (1/mu*Y1_i+X_i-E1_i)';
C = A'*B;
[U,~,V] = svd (C,'econ'); 

PT = U*V';
P = PT';
end
