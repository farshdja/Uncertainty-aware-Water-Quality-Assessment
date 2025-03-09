function [selected_PCs, selected_D, rot_PCs,  rotated_D] = UPCA(data, mu, covMatrices_gmm, pro, id, S)
    % Compute mean vector
    mu_t = zeros(1, size(data, 2));
    for i = 1:size(mu, 1)
        mu_t = mu_t + mu(i, :);
    end
    mu_t = mu_t ./ size(mu, 1);

    % Initialize covariance matrix
    K_tt = zeros(size(data, 2), size(data, 2));

    % Compute K_tt
    for i = 1:size(mu, 1)
        m = data(id{i, 1}, :);  % Extract subset of data
        k = pro(i) * ((1 / length(m)) * (m' * m) + (S .* covMatrices_gmm(:, :, i)) - (mu_t' * mu_t)); 
        K_tt = K_tt + k;
    end

    % Perform PCA using covariance matrix
    coeff = pcacov(K_tt);
    [~, d] = eig(K_tt);
    D = diag(d);

    % Sort eigenvalues and corresponding eigenvectors
    [D, ~] = sort(D, 'descend');
    selected_PCs = coeff(:, 1:size(find(D >= 1)));
    selected_D = D(1:size(find(D >= 1)));
    rot_PCs = rotatefactors(coeff(:, 1:size(find(D >= 1), 1)));
    rotated_D = diag(rot_PCs' * K_tt * rot_PCs);
end

