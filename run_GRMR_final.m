function [X,err]=run_GRMR_final(smp,M,L1,maxIter,gm,rk_sel,I_WB)
 
 
X=I_WB;
[U,S,~]=svd(X);

for iter=1:maxIter
    Uj=U(:,1:rk_sel);
    Pj=Uj*Uj';
    Z=applyFWR(smp,X);
    param_num=norm(Pj*(applyBWR(smp,M-Z)),'fro')^2;
    param_dem=norm(applyFWR(smp,Pj*(applyBWR(smp,M-Z))),'fro')^2;
    param=(param_num)./(param_dem+eps);
    term1=param*applyBWR(smp,(M-Z)) ;
    term2=2*gm*L1*X;
    Y=X+term1-term2; 
    [U,S,V]=svd(Y,'econ');
    S(rk_sel+1:end,rk_sel+1:end)=0;  
    X=U*S*V';
end
