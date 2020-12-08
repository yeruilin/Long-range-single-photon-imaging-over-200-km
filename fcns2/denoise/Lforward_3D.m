function X=Lforward_3D(P)

[m2,n2,p2]=size(P{1});
[m1,n1,p1]=size(P{2});
[m3,n3,p3]=size(P{3});
if (n2~=n1+1||n3~=n1+1)
    error('dimensions are not consistent')
end
if(m1~=m2+1||m3~=m2+1)
    error('dimensions are not consistent')
end
if(p1~=p3+1||p2~=p3+1)
    error('dimensions are not consistent')
end

m=m2+1;
n=n2;
p=p2;

X=zeros(m,n,p);
X(1:m-1,:,:)=P{1};
X(:,1:n-1,:)=X(:,1:n-1,:)+P{2};
X(:,:,1:p-1)=X(:,:,1:p-1)+P{3};
X(2:m,:,:)=X(2:m,:,:)-P{1};
X(:,2:n,:)=X(:,2:n,:)-P{2};
X(:,:,2:p)=X(:,:,2:p)-P{3};



