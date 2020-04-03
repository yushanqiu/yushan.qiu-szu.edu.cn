function [Distance_norm,entropy,s]=entropyss(A)
LA=zeros(size(A));
for i=1:size(A,1)
    for j=1:size(A,2)
        if(A(i,j)~=0)
            LA(i,j)=log(A(i,j));
        end;
    end;
end;
AZ=zscore(LA);

[U,S,V]=svd(AZ);
s=svd(AZ);
plot(s);


whole=0;
for i=1:size(A,1)
    whole=whole+s(i)^2;
end
f=zeros(size(A,1),1);
for i=1:size(A,1)
    f(i)=s(i)^2/whole;
end
middle=0;
for i=1:size(A,1)
    middle=middle+f(i)*log(f(i));
    entropy(i)=(-1)/log(i)*middle;
end


DK=AZ;
PK=zeros(size(A));
DK_norm=zeros(size(A,1),1);
PK_norm=zeros(size(A,1),1);

 for k=1:size(A,1)
    DK=DK-s(k)*U(:,k)*V(:,k)';
    DK_norm(k)=norm(DK,2);%
    for i=1:size(A,1)
        for j=1:size(A,2)
            sign_rand=rand(1);
            sign_rand(sign_rand>0.5)=1;
            sign_rand(sign_rand<=0.5)=-1;
            PK(i,j)=DK(i,j)*sign_rand;%
        end 
    end
     PK_norm(k)=norm(PK,2);%
 end
 
Distance_norm=DK_norm-PK_norm;
plot(Distance_norm);%Ê£Óà·¶Êı