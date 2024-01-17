function K = updateK(H,V,N)
K_complement=zeros(N,N);
L=eye(N)-1/N*ones(N);
for v=1:V
        K_complement = K_complement*0;
        for k=1:V
            if (k==v) 
                continue;
            end
            K_complement = K_complement + L*H{k}'*H{k}*L;                    
        end
        [K{v}] =K_complement;
end
end