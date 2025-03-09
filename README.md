# Uncertainty-aware-Water-Quality-Assessment
This repository introduces an innovative framework for Water Quality Index (WQI) assessment that directly incorporates data uncertainty to enhance the reliability of water quality evaluations. The core function, UWQI_SA.m, performs a comprehensive analysis of water quality data, integrating Gaussian Mixture Model (GMM) fitting, Uncertainty-aware Principal Component Analysis (UPCA), weight computation, WQI calculation, and Sobol Sensitivity Analysis. The primary objectives of this framework are: 1) Examining the impact of varying levels of uncertainty in water quality data on Principal Component Analysis (PCA) outcomes. 2) Quantifying how weighting uncertainty influences the WQI due to input data variability. 3) Identifying, through sensitivity analysis, which water quality parameters significantly contribute to the propagation of uncertainty in the WQI at each monitoring station.


Usage Guidelines
To effectively employ this framework, follow these steps:

1) Data Preprocessing: Utilize the data_preprocessing.m function to prepare your water quality data. This step ensures that the data is cleansed of outliers and standardized to a consistent temporal resolution, facilitating accurate analysis.

2) Uncertainty-aware Analysis: Execute the UWQI_SA.m function with the preprocessed data. This function will conduct the GMM fitting, UPCA, weight computations, WQI calculations, and Sobol Sensitivity Analysis. The outcomes will be stored in the output structure for further examination and interpretation. To utilize and test the proposed approach, ensure that MATLAB is installed on your system. For detailed system requirements and installation guidelines, please refer to the official MATLAB documentation.
