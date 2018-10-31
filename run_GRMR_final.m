function [X,err]=run_GRMR_final(smp,M,L1,maxIter,gm,rk_sel,I_WB)
 
X=I_WB;
[U,~,~]=svd(X);

for iter=1:maxIter
    Uj=U(:,1:rk_sel);
    Pj=Uj*Uj';
    Z=applyFWR(smp,X);
    Q=applyBWR(smp,M-Z);
    W=Pj*Q;
    param_num=norm(W,'fro')^2;
    param_dem=norm(applyFWR(smp,W),'fro')^2;
    param=(param_num)./(param_dem+eps);
    term1=param*Q ;
    term2=2*gm*L1*X;
    Y=X+term1-term2; 
    [U,S,V]=svd(Y,'econ');
    S(rk_sel+1:end,rk_sel+1:end)=0;  
    X=U*S*V';
    if sqrt(norm(Y-X))<3
       break 
    end
end
