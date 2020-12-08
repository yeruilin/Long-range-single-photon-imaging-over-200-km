function [ kNN_Init ] = fcn_kNN_Init( I_tv, T_ML, numK )
% [ kNN_Init ] = fcn_kNN_Init( I_tv, T_ML, numK )
% Use k nearest neighbors to initialize depth regularization

[Lr,Lc] = size(I_tv);

% choose larger dimension (row or column) to standardize position data
dims = [Lr,Lc];
[~, max_dim] = max(dims);
std_dev_dim = std(1:dims(max_dim));

[px_col2, px_row2] = meshgrid(([1:Lc]-(Lc+1)/2)/std_dev_dim,...
    ([1:Lr]-(Lr+1)/2)/std_dev_dim);
[true_col, true_row] = meshgrid(1:Lc,1:Lr);
true_coords = [true_col(:),true_row(:)];

% standardize reflectivity data
px_refl = (I_tv-mean(I_tv(:)))/std(I_tv(:));

Px_data2 = [px_col2(:),px_row2(:),px_refl(:)];
Px_label2 = T_ML(:);

data_train2 = Px_data2(T_ML(:)~=0,:);
label_train2 = Px_label2(T_ML(:)~=0,:);
data_test2 = Px_data2(T_ML(:)==0,:);
true_test = true_coords(T_ML(:)==0,:);

nearest_pixels = knnsearch(data_train2,data_test2,'K',numK+1);

kNN_Init = T_ML;
for ii = 1:length(data_test2)
    kNN_Init(true_test(ii,2),true_test(ii,1)) = ...
        median(label_train2(nearest_pixels(ii,2:end)));
end

end

