function [ Sobol_out ] = Sobol_SA(N, lb, ub, st)

global station
station = st;
% % ****       Sobol' Variance-based Global Sensitivity Analysis        ****
% % ****     This numerical (Monte Carlo-based) integration method      ****
% % ****            is adopted from  Saltelli et al. (2008)             ****
% % ****  where page 164 states "This procedure is the best available   ****
% % **** today for computing indices based purely on model evaluations" ****
% % ************************************************************************
% %% Control Parameters and Specifications
% N = 100;                              % 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
seedNum = [];                  % 2. Seed for random number generator; put [] for autmatic randomization 
 funPath = 'C:\Users\Lenovo\Desktop\software';
 funFile = 'funcWQI_Sobol';            % 4. Model/function file: MATLAB m-file without .m extension
 smplMtd = 'LHS';                     % 5. Sampling Method: RND, LHS, SymLHS, PLHS, SobolSeq, or Halton for generation of star centers; if blank, default is LHS

 bootstrapFlag = 0;                    % 7. Bootstrap Flag: enter 1 to bootstrap, 0 not to bootstrap
 bootstrapSize = 1000;                 % 8. Bootstrap Size: number of sampling iterations with replacement - if bootstrap flag = 0, this will be ignored.
 confdLvl = 0.9;                       % 9. Confidence Level: to report bootstrap-based lower and upper bounds on VARS results; if bootstrapFlag = 0, this line will be ignored.
 numGrp = 3;                           % 10. User-specified number of groups: if blank, VARS-TOOL will find the optimal number; if bootstrapFlag = 0, this line will be ignored.

%% Store Sobol Control Parameters and Specifications
 Sobol_inp.N = N; 
 Sobol_inp.seedNum = seedNum;
 Sobol_inp.funFile = funFile;
 Sobol_inp.funPath = funPath;
 Sobol_inp.funFile = funFile;
 Sobol_inp.smplMtd = smplMtd;
%% Randomization
 if isempty( seedNum ) == false
     rand('state',seedNum); 
 end

%% Generate the base sample from a unit hypercube

numDim = length(lb);  % number of factors
SobolCost = N * (numDim + 2); % total number of function evaluations (model runs)

switch smplMtd
    case 'RND'
        baseSample = rand(N, numDim * 2); %uniformly distributed random numbers
    case 'LHS'
        baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
    case 'SymLHS'
        baseSample = SymLHS( N , numDim * 2,'maxmin', 10 );
    case 'PLHS'
        numSlices = ceil(N/sliceSize);
        if numSlices > 1
            smpSize = numSlices * sliceSize;
            p = PLHSdesign(smpSize, numDim * 2, numSlices, 10, 'maxmin');
            baseSample = p(1:N,:);
            % generates optimal PLHS based on 'maxmin' criterion.
            % by default the number of iteration is 10.
        else
            baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
        end
    case 'Sobolseq'
        p = sobolset(numDim * 2,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
        baseSample = p(1:N,:);
    case 'Halton'
        p = haltonset(numDim * 2,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'RR2'); %PR2 scramling strategy. For no scrambling disable this line.
        baseSample = p(1:N,:);
    otherwise
        baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
end

%% Define sub-sample matrices A, B, and C, and scale them onto the actual factor space
A = baseSample(:, 1 : numDim);
B = baseSample(:, numDim + 1 : end);
for j = 1 : N
    A(j, :) = A(j, :) .* ( ub' - lb') + lb';
    B(j, :) = B(j, :) .* ( ub' - lb') + lb';
end
for i = 1 : numDim
    C{i} = B;
    C{i}(:, i) = A(:, i);
end
%% Run the function/model for all points in matrices A, B, and C and return the function/model response
[yA, yB, yC ] = funcEval(A, B, C, funPath, funFile);

%% Cacluate Sobol's First-Order (FO) and Total-Order (TO) sensitivity indices and rank factors based on TO
[ FO, TO, V, Mu ] = SobolCalc (yA, yB, yC);
rnkFO = factorRanking(FO);
rnkTO = factorRanking(TO);
Sobol_out.FO = FO;
Sobol_out.TO = TO;
Sobol_out.rnkFO = rnkFO;
Sobol_out.rnkTO = rnkTO;
Sobol_out.V = V;
Sobol_out.Mu = Mu;

%% Bootstrap
if bootstrapFlag == 1;
    [ Sobol_out.bootstrap ] = bootstrapSobol (yA, yB, yC, bootstrapSize, confdLvl, rnkFO, rnkTO, numGrp);
end

%% Store Results
% save ('Results_Sobol', 'Sobol_out');
end

%% ************************************************************************
%  ***********                  Sub-functions                   ***********
%  ******                                                            ******
%  **                                                                    **
%% ************************************************************************
function [ FO, TO, Vy, MUy ] = SobolCalc (yA, yB, yC)
[ N,  numDim ] = size(yC);
f0 = mean([ yA; yB ]);  % sample mean - approximation of the mean of response surface
MUy = f0;
Vy = ( [ yA; yB ]' * [ yA; yB ] ) / ( 2 * N ) - f0 ^ 2; % sample variance - approximation of the variance of response surface
for i = 1 : numDim
    FO(i) = ( (( yA' * yC(:, i) ) / N) - mean(yA) * mean(yB) ) / Vy;     % First-Order effect
    TO(i) = 1 - ( ( yB' * yC(:, i) ) / N - f0 ^ 2 ) / Vy; % Total-Order effect
end
end
%% ************************************************************************
function [ yA, yB, yC ] = funcEval (A, B, C, funPath, funFile)
currDir = pwd;
cd(funPath);
[ N , numDim ] = size(A);
SobolCost = N * (numDim + 2); % total number of function evaluations (model runs)

yA = zeros(N, 1);
yB = zeros(N, 1);
yC = zeros(N, numDim);
fprintf('A Sobol experiment started: size of base sample = %g, total number of model runs = %g. \n', N, SobolCost );
for j = 1 : N
    fprintf('Group run (base sample) #%g started. Running model %s %g times...', j, funFile, numDim + 2 );
    tic;
    yA(j, 1) = feval(funFile, A(j, :) );
    yB(j, 1) = feval(funFile, B(j, :) );
    for i = 1 : numDim
        yC(j, i) = feval(funFile, C{i}(j, :) );
    end
    time = toc;
    fprintf(' Group run finished in %g seconds.\n', time);
end
cd (currDir);
end
%% ------------------------------------------------------------------------
function [ bootstrap ] = bootstrapSobol (yA, yB, yC, bootstrapSize, confdLvl, rnkFOBnchmrk, rnkTOBnchmrk, numGrp)

[ N  numDim ] = size(yC);
[randstream1] = RandStream.create('mrg32k3a','NumStreams',1, 'Seed', 1234567);
for k = 1 : bootstrapSize
    rows = ceil ( rand(randstream1, N, 1) * N );
    yAboot = yA(rows);
    yBboot = yB(rows);
    yCboot = yC(rows, :);
    [ FO(k, :), TO(k, :), V(k, :), MU(k, :) ] = SobolCalc (yAboot, yBboot, yCboot);
     rnkFO(k, :) = factorRanking(FO(k, :));
     rnkTO(k, :) = factorRanking(TO(k, :));
end
FO_sorted = sort(FO, 1);
TO_sorted = sort(TO, 1);

bootstrap.FO_low = FO_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.FO_upp = FO_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

bootstrap.TO_low = TO_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.TO_upp = TO_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

for D = 1 : numDim
    bootstrap.rel_FO(1, D) = length ( find( rnkFO(:, D) == rnkFOBnchmrk(1, D) ) ) / bootstrapSize;
    bootstrap.rel_TO(1, D) = length ( find( rnkTO(:, D) == rnkTOBnchmrk(1, D) ) ) / bootstrapSize;
end
% ************ based on group ranking 
rank_benchmark_grp = groupRanking ( rnkTOBnchmrk , numGrp );
for iter = 1 : bootstrapSize
    rnkTO_grp( iter, :) = groupRanking ( rnkTO(iter, :) , numGrp );
end
for D = 1 : numDim
    relGrp(1, D) = length ( find( rnkTO_grp(:, D) == rank_benchmark_grp(1, D) ) ) / bootstrapSize;
end
 
end
%% ************************************************************************
function [ rank ] = factorRanking(SAindices)
[sorted, order] = sort(SAindices, 'descend');
temp = [ order; 1: length(SAindices) ]';
temp2 = sortrows(temp, 1)';
rank = temp2(2, :);
end
%% ************************************************************************
function rank_grp = groupRanking ( rank_indvl, numGrp )
numDim = length (rank_indvl);
grpSize = round ( numDim / numGrp );
grpNum = 1; temp = 0;
for rankNum = 1 : numDim
    if temp == grpSize
        temp = 0;
        grpNum = grpNum + 1;
    end
    temp = temp + 1;
    rank_grp ( 1,  rank_indvl == rankNum  ) = grpNum;
end
end

function [ WQI_out ] = funcWQI_Sobol(x)

FileName   = 'Sub_index.mat';
FolderName = 'C:\Users\Lenovo\Desktop\software';
File       = fullfile(FolderName, FileName);
load(File);   % not: load('File')
global station
WQI_out = sum(x.*Q(station,:))./sum(x);

end