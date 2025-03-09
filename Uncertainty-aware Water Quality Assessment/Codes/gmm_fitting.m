function [pro, mu, covMatrices, id] = gmm_fitting(data, numComponents)
    rng(0); % For reproducibility

    % Set options for GMM fitting
    options = statset('Display', 'final', 'MaxIter', 10000);

    % Fit Gaussian Mixture Model (GMM)
    gmm = fitgmdist(data, numComponents, 'RegularizationValue', 0.01, 'Options', options);

    % Compute outputs
    likelihoods = sum(log(pdf(gmm, data)));
    AICs = gmm.AIC;
    BICs = gmm.BIC;
    pro = gmm.ComponentProportion;
    mu = gmm.mu;
    covMatrices = gmm.Sigma;
    clus = cluster(gmm,data);
    for i = 1:max(clus)
        id{i,1} = find(clus==i);
    end
end

