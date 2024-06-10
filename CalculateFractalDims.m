% CalculateFractalDims

% This script calculates the fractal dimension (i.e., the correlation
% dimension) for each path through the tree of life generated either from
% scripts "BuildMultifractalTree.m" or
% "BuildMultifractalTreeForVisual.m". The script calculates the
% correlation dimension of each path at each iterate.

% Generate a vector of the logarithmically spaced epsilon values
ep_Orig = [logspace(10, 1, 100) / 1.0e+10, logspace(10, 1, 100) /...
    1.0e+20, logspace(10, 1, 100) / 1.0e+30, logspace(10, 1, 100) /...
    1.0e+40, logspace(10, 1, 400) / 1.0e+50, logspace(10, 1, 400) / 1.0e+60];
ep_Orig = ep_Orig(2:end);
for k = 1:size(S, 2)
    for p = 2:size(S, 1)
        
        % For each epsilon value determine the number of pairs of points
        % whose separation distance is less than epsilon
        for i = 1:size(ep_Orig, 2)
            All_Combos = nchoosek(S(1:p, k), 2);
            Differences = All_Combos(:, 1) - All_Combos(:, 2);
            
            % Store the number in C_ep
            C_ep(i) = size(Differences(Differences < ep_Orig(i)), 1);
        end
        
        C_ep = C_ep / p^2; % Divide by the number of points squared
        
        How_Big = size(C_ep(C_ep ~= 0), 2); % Determine how many entries
                                            % in C_ep are not equal to zero
                                            
        C_ep = C_ep(1:How_Big); % Eliminate entries equal to zero.
        
        ep_Altered = ep_Orig(1:How_Big); % Get a truncated version of
                                         % ep_orig to match the truncated
                                         % version of C_ep
                                         
        Lin_Model = fitlm(log(ep_Altered), log(C_ep)); % Fit a linear model
        
        D_F(p, k) = Lin_Model.Coefficients{2, 1}; % The fractal dimension
                                                  % is equal to the slope.
                                                  
        % Get the coefficient of determination for each path.                                          
        R_Squared(p, k) = Lin_Model.Rsquared.Ordinary; 
    end
end

% Get the average R-squared values for the paths at all iterates, as well
% as the standard deviations
Avg_R2_Vector = zeros(1, size(S, 1));
Std_Devs = zeros(1, size(S, 1)); 
for i = 1:size(R_Squared, 1)
    temp = unique(R_Squared(i, :), 'stable');
    Std_Devs(i) = std(temp);
    Avg_R2_Vector(i) = sum(temp) / size(temp, 2);
end

% Look at the fit for a randomly chosen path
Random_Member = randi(size(S, 2)); % Choose a random path
for i = 1:size(ep_Orig, 2)
    All_Combos = nchoosek(S(:, Random_Member), 2);
    Differences = All_Combos(:, 1) - All_Combos(:, 2);
    C_ep_Rand_Member(i) = size(Differences(Differences < ep_Orig(i)), 1);
end
C_ep_Rand_Member = C_ep_Rand_Member / p^2; % Divide by the number of points squared
How_Big = size(C_ep_Rand_Member(C_ep_Rand_Member ~= 0), 2); % eliminate 0's
C_ep_Rand_Member = C_ep_Rand_Member(1:How_Big); % Eliminate entries equal to zero.
ep_Rand_Member = ep_Orig(1:How_Big); % truncate
Lin_Model = fitlm(log(ep_Rand_Member), log(C_ep_Rand_Member)); % Regression
D_F_Rand_Member = Lin_Model.Coefficients{2, 1}; % Get fractal dimension
R_Squared_Rand_Member = Lin_Model.Rsquared.Ordinary; % Get the coefficient
                                                     % of determination


