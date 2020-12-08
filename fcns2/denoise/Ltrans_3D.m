function P=Ltrans_3D(X)

[m,n,p]=size(X);

P{1}=X(1:m-1,:,:)-X(2:m,:,:);
P{2}=X(:,1:n-1,:)-X(:,2:n,:);
P{3}=X(:,:,1:p-1)-X(:,:,2:p);