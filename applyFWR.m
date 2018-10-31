function Y=applyFWR(smp_mtx,X)

% [n1,n2,n3]=size(smp_mtx);
% for xx=1:n1
%     tmp1=squeeze(X(xx,:));
%     for tt=1:n3
%         tmp2=squeeze(smp_mtx(xx,:,tt));
%         Y(xx,tt)=trace(tmp2'*tmp1);
%     end
% end

[n1,n2,n3]=size(smp_mtx);
% Y=zeros(n1,n3);
% for xx=1:n1
%     tmp1=X(xx,:);
%     for tt=1:n3
%         tmp2=smp_mtx(xx,:,tt);
%         Y(xx,tt)=sum(tmp2.*tmp1,2);
%     end
% end

Y=zeros(n1,n3);
for qq=1:n3
    a=X';
    b=squeeze(smp_mtx(:,:,qq));
    Y(:,qq)=sum(b.*a',2);
end
