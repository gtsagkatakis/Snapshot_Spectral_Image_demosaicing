function Y=applyFWR(smp_mtx,X)

[n1,n2,n3]=size(smp_mtx);
for xx=1:n1
    tmp1=squeeze(X(xx,:));
    for tt=1:n3
        tmp2=squeeze(smp_mtx(xx,:,tt));
        Y(xx,tt)=trace(tmp2'*tmp1);
    end
end