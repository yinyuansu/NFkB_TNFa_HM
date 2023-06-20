function [U, S] = pca_self(X)
%   Returns the eigenvectors U, the eigenvalues (on diagonal) in S of input X

% Useful values
[m, n] = size(X);
X_norm = zeros(m, n);
U = zeros(n);
S = zeros(n);

% Normalize before run PCA (mean=0, standard deviation=1) 

mu = mean(X);
X_norm = bsxfun(@minus, X, mu);

sigma = std(X_norm);
X_norm = bsxfun(@rdivide, X_norm, sigma);

Sigma = (X_norm'*X_norm)/m;
[U, S, V] = svd(Sigma);

end
