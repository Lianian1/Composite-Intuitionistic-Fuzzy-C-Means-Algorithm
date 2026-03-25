function M = kmeans_plus_plus_init(X, k)
    % X is the data matrix (Nxd)
    % k is the number of clusters

    [N, ~] = size(X);
    M = zeros(k, size(X, 2)); % Preallocate matrix for cluster centers

    % Randomly select the first cluster center
    M(1, :) = X(randi(N), :);

    for i = 2:k
        D = pdist2(X, M(1:i-1, :), 'squaredeuclidean');
        D = min(D, [], 2);
        probs = D / sum(D);
        cumprobs = cumsum(probs);
        r = rand();
        idx = find(cumprobs >= r, 1, 'first');
        M(i, :) = X(idx, :);
    end
end

% Usage:
% M = kmeans_plus_plus_init(X, k);
