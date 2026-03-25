function j_fun = object_fun(N, d, k, Cluster_elem, landa, M, fuzzy_degree, W, z, beta_z, p, X, v)
    % 特征选择或变换（这里使用 PCA）
    [~, ~, ~, X_transformed, ~] = pca(X);

    C = cov(X_transformed);  % 计算协方差矩阵
    C_inv = inv(C);          % 计算协方差矩阵的逆

    dNK = zeros(N, k);

    for j = 1:k
        diff = X_transformed - repmat(M(j, :), N, 1);
        distance = zeros(N, 1);

        for i = 1:N
            % 计算马氏距离
            distance(i) = sqrt((diff(i, :)')' * C_inv * (diff(i, :)'));
        end

        distance = (1 - exp(-landa .* distance.^2)) .* v;
        WBETA = transpose(z(j, :).^beta_z);
        WBETA(WBETA == inf) = 0;
        dNK(:, j) = distance * WBETA * W(1, j)^p;
    end

    % 引入自适应权重
    adaptive_weight = std(X_transformed);  % 使用标准差作为权重
    dNK = dNK .* adaptive_weight;

    % 计算目标函数值
    j_fun = sum(sum(dNK .* transpose(Cluster_elem.^fuzzy_degree)));
end
