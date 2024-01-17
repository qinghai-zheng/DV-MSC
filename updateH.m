function H_i = updateH(X_i,P_i,Z_i,E1_i,E2_i,Y1_i,Y2_i,K_i,mu,rho,lambda_1,N)


    %D{v} = D{v}./repmat(sqrt(sum(D{v}.^2,1)),size(D{v},1),1);
    I = eye(N,N);
    A = mu*(P_i'*P_i);
    B = lambda_1*(K_i+K_i')-rho*(Z_i+Z_i'-Z_i*Z_i'-I);
    C = -(Y2_i-Y2_i*Z_i'-P_i'*Y1_i-mu*P_i'*(X_i-E1_i)-rho*E2_i+rho*E2_i*Z_i');
    H_i = lyap(A,B,C);
    
end  
