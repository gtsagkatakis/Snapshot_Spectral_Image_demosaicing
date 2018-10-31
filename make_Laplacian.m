function L1 = make_Laplacian(M,sgm2)
    
    W=eye(size(M,1));
    
    for xx=1:size(W,1)
        for yy=1:size(W,1)
            tmp1=M(xx,:);
            tmp2=M(yy,:);
            W(xx,yy)=exp(-(norm(tmp1-tmp2))/(sgm2^2));
        end
    end
    
    DG=sum(W);
    DGmn=DG.^-0.5;
    
    L1=diag(DGmn)*(diag(DG)-W)*diag(DGmn);