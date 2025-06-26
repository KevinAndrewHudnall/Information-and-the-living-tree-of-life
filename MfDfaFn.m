function [Generalized_Hurst_values] = MfDfaFn(S, q_Values, Box_Sizes)

% This function performs Multifractal Detrended Fluctuation Analysis
% (MF-DFA) on a set of time series derived from paths through a
% multifractal tree.
%
% INPUTS:
%   S          - matrix of size [N x M] where each column is a time series
%                representing a path through the tree of life
%   q_Values   - 1D array of moments (q) to analyze (e.g., -5:0.5:5)
%   Box_Sizes  - 1D array of box sizes (scales) used in fluctuation analysis
%
% OUTPUT:
%   Generalized_Hurst_values - [M x length(q_Values)] matrix of H(q) exponents
%                              for each series (column of S) and each q
%
% The function implements the standard MF-DFA algorithm:
%   1. Integrate each detrended time series
%   2. Divide into non-overlapping boxes of varying sizes
%   3. Detrend each box and compute root-mean-square fluctuations
%   4. Aggregate these fluctuations using generalized averaging for each q
%   5. Fit log–log scaling of fluctuation function vs. box size
%
% The slope of the scaling curve gives the generalized Hurst exponent H(q),
% capturing the degree of multifractality in each lineage path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the size of the matrix S
[N, Num_Series] = size(S);

% Initialize a matrix to store the generalized fluctuation function for
% each time series, q, and box size
MF_Fq_Values = zeros(Num_Series, length(Box_Sizes), length(q_Values));

% Loop over each column in S
for j = 1:Num_Series
    % Extract the j-th column
    X = S(:, j);
    
    % Remove the mean from the time series to make it zero-mean
    X = X - mean(X);
    
    % Calculate the cumulative sum
    Y = cumsum(X);
    
    % Loop over each q value
    for q_Idx = 1:length(q_Values)
        q = q_Values(q_Idx);  % Current q value
        
        % Loop over each box size
        for Idx = 1:length(Box_Sizes)
            n = Box_Sizes(Idx);
            Num_Boxes = floor(N / n);  % Number of boxes
            
            % Initialize the fluctuation function for the current box size
            F_nq = zeros(Num_Boxes, 1);
            
            % Loop over each box
            for i = 1:Num_Boxes
                % Define the indices for the current box
                Start_Idx = (i-1) * n + 1;
                End_Idx = i * n;
                
                % Extract the data for the current box
                Box_Data = Y(Start_Idx:End_Idx);
                
                % Fit a linear Trend to the data in the current box
                t = (1:n)';
                p = polyfit(t, Box_Data, 1);  % Linear fit
                Trend = polyval(p, t);
                
                % Calculate the detrended fluctuation for the current box
                F_n = sqrt(mean((Box_Data - Trend).^2));  % RMS of residuals
                
                % Calculate generalized fluctuation function for each q
                if q == 0
                    % For q = 0, use log averaging to avoid singularities
                    F_nq(i) = exp(0.5 * mean(log(F_n.^2)));
                else
                    F_nq(i) = (mean(F_n.^q))^(1/q);
                end
            end
            
            % Avg the generalized fluctuation function over all boxes of size n
            MF_Fq_Values(j, Idx, q_Idx) = mean(F_nq);
        end
    end
end

% Perform log-log fitting to estimate the generalized Hurst exponent 
% for each column of S and q-value
Generalized_Hurst_values = zeros(Num_Series, length(q_Values));

for j = 1:Num_Series
    for q_Idx = 1:length(q_Values)
        log_Box_Sizes = log(Box_Sizes);

        % Extract Fq values for this series and q
        log_Fq = log(squeeze(MF_Fq_Values(j, :, q_Idx)));  
        
        % Fit a linear model to the log-log plot
        p = polyfit(log_Box_Sizes, log_Fq, 1);
        
        % The slope of the log-log plot is the generalized Hurst exponent
        % for the j-th series and q-th moment
        Generalized_Hurst_values(j, q_Idx) = p(1);
    end
end
end