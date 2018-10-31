function Y=applyBWR(smp_mtx,X)

% [n1,n2,n3]=size(smp_mtx);
% Y=zeros(n1,n2);
% for tt=1:n1
%     
%     for qq=1:n3
%         tmp2=X(tt,qq);
%         tmp1=squeeze(smp_mtx(tt,:,qq));
%         Y(tt,:)=Y(tt,:)+(tmp1'*tmp2)';  
%     end
% end

[n1,n2,n3]=size(smp_mtx);
Y=zeros(n1,n2);
% tmp1=zeros(1,n2);
% for tt=1:n1
%     for qq=1:n3
% %         tmp2=X(tt,qq);
%         tmp1=smp_mtx(tt,:,qq);
% %         tmp1=squeeze(smp_mtx(tt,:,qq));
% %         Y(tt,:)=Y(tt,:)+(tmp1'*tmp2)';  
%         Y(tt,:)=Y(tt,:)+(tmp1'*X(tt,qq))';  
%     end
% end

for ww=1:n1
    a=squeeze(smp_mtx(ww,:,:));
    b=squeeze(X(ww,:));
    Y(ww,:)=a*b';
end