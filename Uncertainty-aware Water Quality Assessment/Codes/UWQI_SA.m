
function [output] = UWQI_SA(data, numComponents, s, rand_num, st, N)

% UWQI_SA Uncertainty-Aware Water Quality Index-Sensitivity Analysis
%
% This function performs a comprehensive analysis of water quality data,
% incorporating Gaussian Mixture Model (GMM) fitting, Uncertainty-aware
% Principal Component Analysis (UPCA), weight computation, Water Quality
% Index (WQI) calculation, and Sobol Sensitivity Analysis. The primary
% objectives are:
% 1. Examine the effect of different levels of uncertainty in water quality
%    data on the results of the principal component analysis.
% 2. Quantify the impact of parameter uncertainty on the WQI due to the
%    uncertainty of the input data.
% 3. Identify, through sensitivity analysis, which quality parameters
%    significantly contribute to the propagation of uncertainty in the WQI
%    at each monitoring station.
%
% INPUTS:
%   data          - Spatiotemporal multidimensional preprocessed water
%                   quality data (matrix).
%   numComponents - Number of Gaussian mixtures fitted to the water quality
%                   data, considered as uncertainty ellipses in this study
%                   (integer).
%   s             - Vector determining the uncertainty levels for adjusting
%                   the uncertainty, with values ranging from 0 to 1 (vector).
%   rand_num      - Number of random samples from the empirical distribution
%                   of data obtained using the Locally Adaptive Kernel
%                   Bandwidth Optimization Method (Shimazaki and Shinomoto,
%                   2010) to perform uncertainty analysis via the Monte Carlo
%                   method (integer).
%   st            - Identification Number of station considered for performing the Sobol
%                   sensitivity analysis to determine the parameters affecting
%                   the propagation of uncertainty at each monitoring station
%                   (integer).
%   N             - Number of samples of weights extracted using the Latin
%                   Hypercube method for performing the Sobol sensitivity
%                   analysis (integer).
%
% OUTPUT:
%   output - Structure containing the following fields:
%       .GMM           - Structure with GMM fitting results:
%           .proportion   - Proportions of each Gaussian component.
%           .mu           - Means of each Gaussian component.
%           .covMatrices  - Covariance matrices of each Gaussian component.
%           .Clusters     - Cluster assignments for each data point.
%       .PCA           - Structure with PCA results:
%           .selected_PCs - Selected principal components for each uncertainty level.
%           .Eigenvalues  - Eigenvalues corresponding to the selected principal components.
%           .Varimaxrotated - Varimax rotated principal components.
%           .Varimaxrotatedeigen - Eigenvalues of the varimax rotated components.
%       .Weights       - Structure with weight computation results:
%           .Upperbound_W1 - Upper bounds of weights for method W1.
%           .Lowerbound_W1 - Lower bounds of weights for method W1.
%           .Upperbound_W2 - Upper bounds of weights for method W2.
%           .Lowerbound_W2 - Lower bounds of weights for method W2.
%       .UWQI_95ppu    - Structure with 95% predictive uncertainty bounds for WQI:
%           .Upperbound_WQI1 - Upper bounds of WQI for method WQI1.
%           .Lowerbound_WQI1 - Lower bounds of WQI for method WQI1.
%           .Upperbound_WQI2 - Upper bounds of WQI for method WQI2.
%           .Lowerbound_WQI2 - Lower bounds of WQI for method WQI2.
%       .Sobol         - Structure with Sobol sensitivity analysis results.
%
% REFERENCES:
%   Gortler, J., Spinner, T., Streeb, D., Weiskopf, D., Deussen, O., 2020. Uncertainty-Aware Principal Component Analysis. IEEE Trans. Visual. Comput. Graphics 26, 822–831.
%   Shimazaki, H., Shinomoto, S., 2010. Kernel bandwidth optimization in spike rate estimation. J Comput Neurosci 29, 171–182.
%
% EXAMPLE USAGE:
%   [output] = UWQI_SA(data, 3, 0:0.01:1, 1000, 10, 500);
%



%% Gaussian Mixture Model Fitting

[pro, mu_gmm, covMatrices_gmm, id] = gmm_fitting(data, numComponents);
GMM.proportion=pro;
GMM.mu=mu_gmm;
GMM.covMatrices=covMatrices_gmm;
GMM.Clusters=id;
output.GMM = GMM;

%% Uncertainty-aware PCA and Weight Calculation Across Different Levels of Uncertainty

for i = 1:length(s)
%
[selected_PCs{i}, selected_D{i}, rot_PCs{i},  rotated_D{i}] = UPCA(data, mu_gmm, covMatrices_gmm, pro, id, s(i));
%
[W1(:,i), W2(:,i)] = computeWeights(selected_PCs{i}, selected_D{i}, rot_PCs{i}, rotated_D{i});

end
PCA.selected_PCs = selected_PCs;
PCA.Eigenvalues = selected_D;
PCA.Varimaxrotated = rot_PCs;
PCA.Varimaxrotatedeigen = rotated_D;
output.PCA = PCA;
%% 
[EMP_rand_W1, ub_W1, lb_W1] = rand_weights_ssvkernel(W1, rand_num);
[EMP_rand_W2, ub_W2, lb_W2] = rand_weights_ssvkernel(W2, rand_num);
Weights.Upperbound_W1 = ub_W1;
Weights.Lowerbound_W1 = lb_W1;
Weights.Upperbound_W2 = ub_W2;
Weights.Lowerbound_W2 = lb_W2; 
output.Weights = Weights;
%% WQI Calculation

for i = 1:rand_num
 WQI1(i,:) = WQI(EMP_rand_W1(:,i));
 WQI2(i,:) = WQI(EMP_rand_W2(:,i));
end

for i = 1:size(WQI1,2)

[lower_bound_WQI1(i), upper_bound_WQI1(i)] = computeUncertainty(WQI1(:,i));
[lower_bound_WQI2(i), upper_bound_WQI2(i)] = computeUncertainty(WQI2(:,i));

end
UWQI_95ppu.Upperbound_WQI1 = upper_bound_WQI1;
UWQI_95ppu.Lowerbound_WQI1 = lower_bound_WQI1;
UWQI_95ppu.Upperbound_WQI2 = upper_bound_WQI2;
UWQI_95ppu.Lowerbound_WQI2 = lower_bound_WQI2;

output.UWQI_95ppu = UWQI_95ppu;

function [lower_bound, upper_bound] = computeUncertainty(data)
    % Compute 95% PPU bounds
    lower_bound = prctile(data, 2.5);
    upper_bound = prctile(data, 97.5);
end

%% Sobol Sensitivity Analysis
[ Sobol_out ] = Sobol_SA(N, lb_W1', ub_W1', st);
output.Sobol=Sobol_out;

end
