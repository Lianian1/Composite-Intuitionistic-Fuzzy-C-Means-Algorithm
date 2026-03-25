% Specify the path to your dataset
data_path = 'C:\Users\amrian\Desktop\archive(1)\pathbased.txt'; % Change this to your actual file path

% Load data
data = load(data_path);

% Assuming that the data is in two columns for x and y
X = data(:, 1:2);

% Number of clusters
num_clusters = 3; % Change this to the number of clusters you want

% Apply Fuzzy C-Means clustering
[center, U] = fcm(X, num_clusters);

% Find the cluster membership for each data point
% ... After your clustering code
[~, membership] = max(U); % FCM membership matrix transformed into a cluster label vector

% Ensure 'membership' is a column vector
membership = membership(:); 

% If X has been transformed or reduced (e.g., PCA), make sure to use the transformed data
% that corresponds to the clustering solution when calling 'evalclusters'.

% Now, calculate the Davies-Bouldin Index
db_index = evalclusters(X, membership, 'DaviesBouldin');

% Print the Davies-Bouldin Index
fprintf('Davies-Bouldin Index: %f\n', db_index.CriterionValues);

% Plot the clustered data with ellipses
figure;
hold on;
colors = ['r', 'g', 'b']; % Colors for the clusters
for i = 1:num_clusters
    % Select data points based on the cluster membership
    cluster_points = X(membership == i, :);
    scatter(cluster_points(:, 1), cluster_points(:, 2), colors(i)); % Plot points for each cluster
    
    % Calculate and plot the ellipse for each cluster
    plot_ellipse(cluster_points(:, 1), cluster_points(:, 2), colors(i));
end

% Plot the cluster centers
scatter(center(:, 1), center(:, 2), 'kx'); % 'kx' denotes black cross for the centers

% Enhance the plot
title('FCM Clustering Results with Ellipses');
xlabel('\fontname{Times New Roman}x \it\fontname{Times New Roman}x\rm')
ylabel('\fontname{Times New Roman}x')
legend(arrayfun(@(x) sprintf('Cluster %d', x), 1:num_clusters, 'UniformOutput', false));
grid on;
hold off;

% Calculate clustering evaluation metrics
silhouette_val = silhouette(X, membership);
db_index = evalclusters(X, membership, 'DaviesBouldin');
ch_index = evalclusters(X, membership, 'CalinskiHarabasz');

% Print the evaluation metrics
fprintf('Silhouette Coefficient: %f\n', mean(silhouette_val));
fprintf('Davies-Bouldin Index: %f\n', db_index.CriterionValues);
fprintf('Calinski-Harabasz Index: %f\n', ch_index.CriterionValues);



% Parameters for the data
n_samples = 900000;
n_clusters = 5;
sub_clusters = 4; % Number of sub-clusters within each main cluster
std_dev = 0.5; % Standard deviation of clusters

% Initialize data array
data = [];

% Define a spread factor to ensure clusters are spread out in space
spread_factor = 8;

% Create data for each cluster
for k = 1:n_clusters
    % Generate random cluster center for each cluster
    center = rand(1, 3) * spread_factor * k;
    
    % Generate points around the center
    cluster_data = bsxfun(@plus, randn(n_samples/n_clusters, 3) * std_dev, center);
    
    % Append to the data array
    data = [data; cluster_data];
end

% Apply Fuzzy C-Means clustering
[centers, U] = fcm(data, n_clusters);

% Determine the maximum membership value for each data point
[~, membership] = max(U);

% Initialize figure
figure;
hold on;

% Define colors for each main cluster
main_colors = lines(n_clusters);

% Define colors for sub-clusters within each main cluster
sub_colors = jet(sub_clusters);

% Plot each main cluster with its sub-clusters
for i = 1:n_clusters
    % Find data points for each main cluster
    idx = membership == i;
    cluster_data = data(idx, :);
    
    % Perform sub-clustering within this main cluster
    [sub_idx, ~] = kmeans(cluster_data, sub_clusters, 'Replicates', 3);
    
    % Plot data points for each sub-cluster
    for j = 1:sub_clusters
        sub_cluster_idx = sub_idx == j;
        scatter3(cluster_data(sub_cluster_idx, 1), cluster_data(sub_cluster_idx, 2), cluster_data(sub_cluster_idx, 3), ...
                 [], sub_colors(j, :), 'filled');
    end
end

% Plot the main cluster centers
scatter3(centers(:,1), centers(:,2), centers(:,3), 100, 'kx', 'LineWidth', 2);

% Release hold on figure
hold off;

% Add title and labels
title('3D Fuzzy C-Means Clustering with Sub-Clusters');
xlabel('X\rm axis');
ylabel('Y axis');
zlabel('Z axis');

% Set grid and axis properties
grid on;
axis equal; % Set equal scaling on all axes to maintain aspect ratio

% Set the 3D view angle
view(3); % Sets the default three-dimensional view with azimuth -37.5 and elevation 30

% 初始化用于保存隶属度的矩阵
membership_degrees = zeros(sub_clusters, n_clusters);

for i = 1:n_clusters
    cluster_data = data(membership == i, :);
    [sub_idx, sub_centers] = kmeans(cluster_data, sub_clusters, 'Replicates', 10);
    
    % 计算协方差矩阵及其逆
    cluster_cov = cov(cluster_data);
    cluster_cov_inv = inv(cluster_cov);
    
    for j = 1:sub_clusters
        sub_cluster_data = cluster_data(sub_idx == j, :);
        
        % 计算马氏距离
        diffs = bsxfun(@minus, sub_cluster_data, sub_centers(j, :));
        mahal_dist = sqrt(sum((diffs * cluster_cov_inv) .* diffs, 2));
        
        % 使用高斯函数计算隶属度
        sigma = std(mahal_dist); % 使用马氏距离的标准差
        gaussian_membership = exp(-(mahal_dist.^2) / (2 * sigma^2));
        
        % 计算平均隶属度
        membership_degrees(j, i) = mean(gaussian_membership);
    end
end

% 标准化隶属度
membership_degrees = bsxfun(@rdivide, membership_degrees, sum(membership_degrees));

% Now membership_degrees is a 4x5 matrix containing the membership degrees
% for the sub-clusters within each of the 5 main clusters

% Display the membership degrees matrix
disp('Membership Degrees Matrix:');
disp(membership_degrees);

% ... [之前的代码] ...

% 根据隶属度计算非隶属度
non_membership_degrees = 1 - membership_degrees;
error_margin=0.05;
% 为每个隶属度值确定5%误差的区间值毕达哥拉斯模糊数
lower_bounds = max(0, membership_degrees - error_margin); % 确保下限不小于0
upper_bounds = min(1, membership_degrees + error_margin); % 确保上限不大于1

% 确定理想最佳解和最差解
ideal_best = max(upper_bounds, [], 1);
ideal_worst = min(lower_bounds, [], 1);

% 初始化距离矩阵
distance_to_best = zeros(1, n_clusters);
distance_to_worst = zeros(1, n_clusters);

% 计算每个子簇到理想最佳解和最差解的距离
for i = 1:n_clusters
    distance_to_best(i) = norm(lower_bounds(:, i) - ideal_best(i), 2) + norm(upper_bounds(:, i) - ideal_best(i), 2);
    distance_to_worst(i) = norm(lower_bounds(:, i) - ideal_worst(i), 2) + norm(upper_bounds(:, i) - ideal_worst(i), 2);
end

% 计算相似度得分
similarity_score = distance_to_worst ./ (distance_to_best + distance_to_worst);

% 将相似度得分转化为排序得分
[~, rank_order] = sort(similarity_score, 'descend');

% 输出特征属性的排序
disp('Feature Attribute Ranking:');
disp(rank_order);


% ... [之前的代码] ...

% 假设的判断矩阵（5x5），对五个主簇进行比较
A = [1, 1/3, 5, 7, 1/5; 
     3, 1, 6, 1/2, 3; 
     1/5, 1/6, 1, 1/3, 2; 
     1/7, 2, 3, 1, 1/7;
     5, 1/3, 1/2, 7, 1];

% 计算权重向量 w（AHP方法）
[V, D] = eig(A);
[max_eigenvalue, max_index] = max(diag(D));
w = abs(V(:,max_index)) / sum(abs(V(:,max_index))); % 使用绝对值确保权重为正数

% 一致性检验
CI = (max_eigenvalue - 5) / (5 - 1);
RI = [0, 0, 0.52, 0.89, 1.11, 1.25]; % 对应于n=5的随机一致性指数
CR = CI / RI(5);
fprintf('一致性比率 CR = %f\n', CR);
if CR < 0.10
    fprintf('判断矩阵具有满意的一致性。\n');
else
    fprintf('判断矩阵的一致性不满意，需要重新评估。\n');
end

% 显示计算出的权重
disp('Calculated Weights for Each Main Cluster (Attribute):');
disp(w);

% 假设 membership_degrees 是一个5x4的隶属度矩阵
% membership_degrees = [ ... ]; % 这应该是之前步骤的输出

% ... [之前的代码] ...

% AHP权重，假设已经计算出来，且其维度为5x1
w = [0.2074; 0.2875; 0.0952; 0.1271; 0.2827]; % 示例权重

% 使用AHP权重调整隶属度得分
% 这里需要将权重w转换为1x5的行向量
w_row = w';

% 现在我们将权重应用于membership_degrees的每一行
% 为此我们需要扩展w_row以匹配membership_degrees的列数
weighted_scores = membership_degrees .* repmat(w_row, size(membership_degrees, 1), 1);

% 计算加权得分的总和，沿着列方向
weighted_sum = sum(weighted_scores, 1);

% 按加权得分排序子簇
[~, sub_cluster_ranking] = sort(weighted_sum, 'descend');

% 输出子簇的排序
disp('Sub-cluster Ranking Based on Weighted Scores:');
disp(sub_cluster_ranking);
