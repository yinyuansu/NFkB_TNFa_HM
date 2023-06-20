function [GroupIndex_pos, ix, B] = kMeansTest(X, k)
% kMeansTest using K-means clustering to partition Matrix X into k clusters,
% return the group index of X by GroupIndex_pos, 
% index order (ix) and number(B) from small to big

N = zeros(k,1);
mainGroup_pos = 1;
[GroupIndex_pos, C_pos, sumD_pos, D_pos] = kmeans(X, k, 'MaxIter', 10000,...
    'Replicates', 10, 'Distance', 'sqeuclidean');


%fprintf('Kmeans_pos group ===\n');
for i = 1:k
    %fprintf('%d: %d (%0.2f%%)\n',i, sum(GroupIndex_pos==i), (sum(GroupIndex_pos==i)/size(X,1))*100);  
    N(i) = sum(GroupIndex_pos==i);
end
[B, ix] = sort(N);

%fprintf('MainGroup_pos: %d\n', ix(k));


%{
fprintf('Kmeans_ori group ===\n');
for i = 1:k
    fprintf('%d: %d (%0.2f%%)\n',i, sum(idx_ori==i), (sum(idx_ori==i)/size(X,1))*100);
end
%}
end