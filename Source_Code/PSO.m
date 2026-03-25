function positions = PSO(X, k, num_particles, lb, ub)
    % PSO algorithm for cluster center optimization

    % PSO 参数
    w = 0.5;
    c1 = 1.5;
    c2 = 1.5;
max_iter_pso = 50;  % 或者你想要的迭代次数
lb = zeros(1, d);  % 适当的维度
ub = ones(1, d);   % 适当的维度
positions = PSO(X, k, num_particles, lb, ub);

    % PSO 初始化
    positions = lb + (ub - lb) .* rand(num_particles, size(lb, 2));
    velocities = rand(size(positions));

    % 找到最优解
    [~, best_idx] = min(objective_function(X, positions, k));

    % 记录每个粒子的最佳位置
    personal_best_positions = positions;
    personal_best_values = objective_function(X, personal_best_positions, k);

    % 记录全局最佳位置
    global_best_position = positions(best_idx, :);
    global_best_value = personal_best_values(best_idx);

    % PSO 主循环
    for iter = 1:max_iter_pso
        % 更新速度和位置
        r1 = rand(size(positions));
        r2 = rand(size(positions));
        velocities = w * velocities + c1 * r1 .* (personal_best_positions - positions) + c2 * r2 .* (repmat(global_best_position, num_particles, 1) - positions);
        positions = positions + velocities;

        % 限制位置在搜索空间内
        positions = max(positions, lb);
        positions = min(positions, ub);

        % 更新最优解
        values = objective_function(X, positions, k);
        improved_idx = values < personal_best_values;
        personal_best_positions(improved_idx, :) = positions(improved_idx, :);
        personal_best_values(improved_idx) = values(improved_idx);

        % 更新全局最优解
        [~, best_idx] = min(personal_best_values);
        global_best_position = personal_best_positions(best_idx, :);
        global_best_value = personal_best_values(best_idx);
    end

    % 最终结果为全局最优解
    positions = global_best_position;
end

function value = objective_function(X, positions, k)
    % 目标函数，用于PSO

    % 计算每个聚类中心到数据点的距离
    distances = pdist2(X, positions);

    % 计算每个数据点到最近聚类中心的距离的平方和
    [~, min_indices] = min(distances, [], 2);
    value = sum((X - positions(min_indices, :)).^2, 'all');
end