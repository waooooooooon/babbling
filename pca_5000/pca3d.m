X=[1,2,3,7,5,6,8;4,5,6,5,6,3,4;1,7,6,9,5,3,7;3,4,5,6,7,8,3;3,4,2,6,7,4,2];
%データ数5,評価項目７
B = X - repmat(mean(X), size(X,1), 1); %平均を引く（縦方向が同じ平均値）
C = cov(B);
[V,D] = eig(C);
[s, index] = sort(diag(D), 'descend');
coeff = V(:,index);
latent = s;
score = B * coeff;

%score=X*coeff;
scatter3(score(:,1),score(:,2),score(:,3));
