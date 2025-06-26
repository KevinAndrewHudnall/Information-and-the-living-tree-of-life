function [D_F, R_Squared, Avg_R2_Vector, Std_Devs, Random_Member, C_ep_Rand_Member,...
    ep_Rand_Member, D_F_Rand_Member, R_Squared_Rand_Member] = CalculateFractalDimsFn(S)
% This script calculates the fractal dimension (i.e., the correlation
% dimension) for each path through the tree of life generated from the
% script "BuildMultifractalTree.m". The script calculates the
% correlation dimension of each path at each iterate.
%
% INPUT:
%   S - scale matrix where each column represents a path through the tree,
%       and each row a generation (from BuildMultifractalTreeFn)
%
% OUTPUTS:
%   D_F                  - matrix of estimated fractal dimensions per path
%   R_Squared            - corresponding R² values for linear fit accuracy
%   Avg_R2_Vector        - average R² per generation
%   Std_Devs             - standard deviation of R² across paths per generation
%   Random_Member        - index of a randomly chosen path
%   C_ep_Rand_Member     - epsilon-neighborhood correlation counts (random path)
%   ep_Rand_Member       - corresponding epsilon values
%   D_F_Rand_Member      - estimated fractal dimension of the random path
%   R_Squared_Rand_Member- R² value for the random path
%
% This function estimates the correlation dimension D₂ using log–log scaling
% of the epsilon-neighborhood function C(ε) across all time series paths S.
% It performs regression for each generation and each path to extract
% the scaling exponent and fit quality (R²). A random path is also selected
% for detailed analysis and diagnostic output.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get minimum spacing between leaves
% Compute the minimum spacing between leaves for each column
temp = min(abs(diff(S)));
% Compute the overall minimum spacing
Min_Spacing = min(temp);
% Compute the adjusted minimum spacing
Min_Spacing = abs(floor(log10(Min_Spacing))) + 2;
 
% Generate a vector of the logarithmically spaced epsilon values
ep_Orig = [logspace(10, 1, 100) / 1.0e+10, logspace(10, 1, 100) /...
    1.0e+20, logspace(10, 1, 100) / 10^Min_Spacing];

D_F = zeros(size(S));
ep_Orig = ep_Orig(2:end);
R_Squared = zeros(size(S));
Std_Devs = zeros(1, size(S, 1));
Avg_R2_Vector = zeros(1, size(S, 1));
parfor k = 1:size(S, 2)
    D_F_Col = zeros(size(S, 1), 1);  % Temporary storage for this column
    R_Squared_Col = zeros(size(S, 1), 1);  % Temporary storage for this column

    for p = 3:size(S, 1) % Start at three since first two do not exist
        
        % For each epsilon value determine the number of pairs of points
        % whose separation distance is less than epsilon
        S_Chunk = S(1:p, k);
        All_Combos = nchoosek(S_Chunk, 2);
        Differences = All_Combos(:, 1) - All_Combos(:, 2);

        % Vectorized computation for C_ep
        C_ep = arrayfun(@(ep) sum(Differences < ep), ep_Orig);

        C_ep = C_ep(C_ep ~= 0); % Eliminate entries equal to zero.
        C_ep = C_ep / p^2; % Divide by the number of points squared
        ep_Altered = ep_Orig(C_ep ~= 0); % Get a truncated version of
                                         % ep_orig to match the truncated
                                         % version of C_ep
                                         
        Lin_Model = fitlm(log(ep_Altered), log(C_ep)); % Fit a linear model
        D_F_Col(p) = Lin_Model.Coefficients{2, 1}; % The fractal dimension
                                                  % is equal to the slope.
        % Get the coefficient of determination for each path.                                          
        R_Squared_Col(p) = Lin_Model.Rsquared.Ordinary;
    end

    % Store results for this column
    D_F(:, k) = D_F_Col;
    R_Squared(:, k) = R_Squared_Col;
end

% Get averages and standard deviations
for k = 1:size(S, 1)
    % Get the average R-squared values for the paths at all iterates,
    % as well as the standard deviations
    Std_Devs(k) = std(unique(R_Squared(k, :)));
    Avg_R2_Vector(k) = mean(unique(R_Squared(k, :)));
end


% Look at the fit for a randomly chosen path
Random_Member = randi(size(S, 2)); % Choose a random path
All_Combos = nchoosek(S(:, Random_Member), 2);
Differences = All_Combos(:, 1) - All_Combos(:, 2);

% Vectorized computation for C_ep_Rand_Member
C_ep_Rand_Member = arrayfun(@(ep) sum(Differences < ep), ep_Orig);

C_ep_Rand_Member = C_ep_Rand_Member / size(S, 1)^2; % Divide by the number of points squared
C_ep_Rand_Member = C_ep_Rand_Member(C_ep_Rand_Member ~= 0); % Eliminate entries equal to zero.
ep_Rand_Member = ep_Orig(C_ep_Rand_Member ~= 0); % truncate
Lin_Model = fitlm(log(ep_Rand_Member), log(C_ep_Rand_Member)); % Regression
D_F_Rand_Member = Lin_Model.Coefficients{2, 1}; % Get fractal dimension
R_Squared_Rand_Member = Lin_Model.Rsquared.Ordinary; % Get the coefficient
                                                     % of determination
end        
