function [W1, W2] = computeWeights(selected_PCs, selected_D, rot_PCs, rotated_D)


% Initialize matrices
gama_L = zeros(size(selected_PCs, 1), length(selected_D));
sigma_gama_L = zeros(size(selected_PCs, 1), 1);

% Compute gama_L and sigma_gama_L
for k = 1:size(selected_PCs, 1)
    for j = 1:length(selected_D)
        gama_L(k, j) = selected_D(j) * abs(selected_PCs(k, j));
    end
    sigma_gama_L(k) = sum(gama_L(k, :));
end

% Compute weight vector W and normalize to obtain W1
W = sigma_gama_L / sum(selected_D);
W1 = W ./ sum(W);

% Find indices of maximum absolute values in each row of rot_PCs
adrs = zeros(size(selected_PCs, 1), 1);
for i = 1:size(selected_PCs, 1)
    [~, adrs(i)] = max(abs(rot_PCs(i, :)));
end

% Initialize W2
W2 = zeros(size(selected_PCs, 1), 1);

% Compute W2 based on relative rotated loadings and relative eigenvalues
for i = 1:length(selected_D)
    ADRS = find(adrs == i);
    if ~isempty(ADRS)
        Rot_load = abs(rot_PCs(ADRS, i));  % Rotated loading
        x = Rot_load ./ sum(Rot_load);     % Relative rotated loading
        rel_Rel_load(ADRS) = x;            % Store relative rotated loading

        y = rotated_D ./ sum(rotated_D);   % Relative eigenvalues
        Rel_eig = y;

        W2(ADRS) = Rel_eig(i) * x;         % Compute W2
    end
end
end

