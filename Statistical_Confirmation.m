%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Statistical_Confirmation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script confirms that reported crosspath quantities
% (e.g., percent divergent, convergent, imaginary, real, unresolved) are
% statistically representative of repeated system simulations.
%
% NOTE: This script is computationally intensive. It performs multiple
% full system simulations and should be run on a machine with parallel
% processing capabilities and adequate RAM. Lowering ITERATIONS will reduce
% computation expense.
%
% This script:
%   - Repeats the entire multifractal tree generation and analysis pipeline
%     for a number of independent runs (e.g., 10)
%   - Samples n random leaf paths per run and evaluates nchoosek(n, 2)
%     pairwise comparisons
%   - For each comparison, computes:
%       * Most recent common ancestor (MRCA)
%       * Crosspath entropy (H), mutual information (I), and distance (d)
%       * Solution to the dilation equation and its classification
%   - Tallies the proportion of:
%       * Coherent (positive d) vs. divergent (negative d) comparisons
%       * Real vs. imaginary solutions
%       * Convergent vs. unresolved or indeterminant cases
%   - Aggregates these statistics across runs and reports mean + std dev
%
% The resulting statistics support the representative validity of values
% reported in Figures 10–12 of the manuscript.

clear;

% Parameters used in the manuscript
ITERATIONS = 30;
MaxOffspring = 3;
MaxGens = 2;

% Number of systems to simulate
Num_Runs = 10;

% Vectors for storing results
Stat_Percent_Coherent = zeros(1, Num_Runs);
Stat_Percent_Divergent = zeros(1, Num_Runs);
Stat_Percent_Real = zeros(1, Num_Runs);
Stat_Percent_Imag = zeros(1, Num_Runs);
Stat_Percent_Convergent = zeros(1, Num_Runs);
Stat_Percent_Not_Resolved = zeros(1, Num_Runs);
Stat_Percent_Indeterminant = zeros(1, Num_Runs);

for u = 1:Num_Runs
    % Build the random tree and get the scale matrix S, the offspring
    % matrix P, and the entropy matrix H.
    [S, P, H] = BuildMultifractalTreeFn(MaxOffspring, MaxGens, ITERATIONS);

    % Randomly select indices for the columns we care about
    n = 10000;  % Number of columns to sample, keeping in mind this will be
                % nchoosek(n, 2) samples.
    S_Final = S(end, :);
    Rand_Locs_In_S = randsample(length(S_Final), n);  % Randomly sample indices

    % Initialize the output matrices with zeros
    D_F = zeros(size(S));

    % Get the minimum spacing between leaves
    temp = min(abs(diff(S)));
    Min_Spacing = min(temp);
    Min_Spacing = abs(floor(log10(Min_Spacing))) + 2;

    % Generate a vector of logarithmically spaced epsilon values
    ep_Orig = [logspace(10, 1, 100) / 1.0e+10, logspace(10, 1, 100) / ...
        1.0e+20, logspace(10, 1, 100) / 10^Min_Spacing];
    ep_Orig = ep_Orig(2:end);

    % Preallocate a cell array to store temporary results for each selected column
    D_F_Temp = cell(1, length(Rand_Locs_In_S));

    % Compute only for the selected columns
    parfor idx = 1:length(Rand_Locs_In_S)
        k = Rand_Locs_In_S(idx);  % Get the actual column index
        D_F_Col = zeros(size(S, 1), 1);  % Temporary storage for this column

        for p = 3:size(S, 1)  % Start at three since first two do not exist
            S_Chunk = S(1:p, k);
            All_Combos = nchoosek(S_Chunk, 2);
            Differences = All_Combos(:, 1) - All_Combos(:, 2);

            % Vectorized computation for C_ep
            C_ep = arrayfun(@(ep) sum(Differences < ep), ep_Orig);
            C_ep = C_ep(C_ep ~= 0);  % Eliminate entries equal to zero
            C_ep = C_ep / p^2;  % Divide by the number of points squared
            ep_Altered = ep_Orig(C_ep ~= 0);  % Truncate ep_Orig to match C_ep

            % Fit a linear model to the log-log data
            Lin_Model = fitlm(log(ep_Altered), log(C_ep));
            D_F_Col(p) = Lin_Model.Coefficients{2, 1};  % Fractal dimension (slope)
        end

        % Store the result in the temporary cell array
        D_F_Temp{idx} = D_F_Col;
    end

    % Collect results from the temporary cell array into D_F
    for idx = 1:length(Rand_Locs_In_S)
        k = Rand_Locs_In_S(idx);
        D_F(:, k) = D_F_Temp{idx};
    end

    % Set a convergence tolerance for the crosspath fractal dimension
    Tolerance = 0.1;
    % Obtain a matrix that gives D_F for those within tolerance
    Conv_Tol = abs(D_F(end, :) - D_F);
    % Logical matrix indicating where Conv_Tol is within the tolerance
    Within_Tolerance = Conv_Tol < (Tolerance * D_F(end, :));
    % Initialize temp with zeros
    temp = zeros(size(Conv_Tol));
    % Assign values from D_F to temp where Within_Tolerance is true
    temp(Within_Tolerance) = D_F(Within_Tolerance);
    for i = size(temp, 1):-1:1
        for j = 1:size(temp, 2)
            if (temp(i, j) == 0)
                temp(1:i, j) = 0;
            end
        end
    end
    Conv_Tol = temp;
    clear temp;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the information quantities for nchoosek(n, 2) randomly    %
    % sampled comparisons and the dilation equations for all convergent   %
    % comparisons                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rand_Combos_LoKs = nchoosek(1:size(Rand_Locs_In_S, 1), 2);
    Rand_Combos_LoCs = [ones(size(Rand_Combos_LoKs, 1), 1) * size(S, 1), ...
        Rand_Locs_In_S(Rand_Combos_LoKs(:, 1)), ones(size(Rand_Combos_LoKs, 1), 1) * ...
        size(S, 1), Rand_Locs_In_S(Rand_Combos_LoKs(:, 2))];

    % Get the parents
    Parent_1 = arrayfun(@(row, col) S(row - 1, col), ...
        Rand_Combos_LoCs(:, 1), Rand_Combos_LoCs(:, 2));
    Parent_2 = arrayfun(@(row, col) S(row - 1, col), ...
        Rand_Combos_LoCs(:, 3), Rand_Combos_LoCs(:, 4));

    % Initialize variables
    MRCA_Locs = zeros(size(Rand_Combos_LoCs, 1), 2);
    MRCA_Values = zeros(size(Rand_Combos_LoCs, 1), 1);
    MRCA_Locs_Reduced = zeros(size(Rand_Combos_LoCs, 1), 2);
    Crosspath_H = zeros(size(Rand_Combos_LoCs, 1), 1);
    Crosspath_I = zeros(size(Rand_Combos_LoCs, 1), 1);
    Crosspath_d = zeros(size(Rand_Combos_LoCs, 1), 1);
    Crosspath_D_F = zeros(size(Rand_Combos_LoCs, 1), 1);
    Dil_Eq = zeros(size(Rand_Combos_LoCs, 1), 1);
    Not_Resolved_Count = 0;
    Indeterminant_Count = 0;

    for ch = 1:size(Rand_Combos_LoCs, 1)
        if Parent_1(ch) == Parent_2(ch)
            Ancestors_Loc = [size(S, 1) - 2, Rand_Combos_LoCs(ch, 2)];
        else
            for i = (size(S, 1) - 1):-1:1
                if S(i, Rand_Combos_LoCs(ch, 2)) == S(i, Rand_Combos_LoCs(ch, 4))
                    Ancestors_Loc = [i, Rand_Combos_LoCs(ch, 2)];
                    break;
                end
            end
        end

        MRCA_Locs(ch, :) = Ancestors_Loc;
        MRCA_Values(ch) = S(MRCA_Locs(ch, 1), MRCA_Locs(ch, 2));

        Crosspath_H(ch) = ((2 * Parent_1(ch) * Parent_2(ch)) / MRCA_Values(ch)^2) * log(MRCA_Values(ch));
        Crosspath_I(ch) = (Parent_1(ch) * Parent_2(ch) / MRCA_Values(ch)^2) * log((Parent_1(ch) * Parent_2(ch)) / MRCA_Values(ch)^2);
        Crosspath_d(ch) = Crosspath_H(ch) - Crosspath_I(ch);

        if Conv_Tol(MRCA_Locs(ch, 1), MRCA_Locs(ch, 2)) ~= 0
            MRCA_Locs_Reduced(ch, :) = MRCA_Locs(ch, :);
            if MRCA_Locs_Reduced(ch, 1) == 1 || MRCA_Locs_Reduced(ch, 1) == 2
                Not_Resolved_Count = Not_Resolved_Count + 1;
            elseif abs(log(S(size(S, 1) - 1, Rand_Combos_LoCs(ch, 2)))) == abs(log(S(size(S, 1) - 1, Rand_Combos_LoCs(ch, 4))))
                Indeterminant_Count = Indeterminant_Count + 1;
            elseif MRCA_Locs_Reduced(ch, 1) ~= 0 && MRCA_Locs_Reduced(ch, 2) ~= 0
                Crosspath_D_F(ch) = D_F(MRCA_Locs_Reduced(ch, 1), MRCA_Locs_Reduced(ch, 2));
                if abs(log(S(size(S, 1) - 1, Rand_Combos_LoCs(ch, 1)))) > abs(log(S(size(S, 1) - 1, Rand_Combos_LoCs(ch, 2))))
                    Dil_Eq(ch) = (1 + (Crosspath_I(ch) / Crosspath_d(ch))^(Crosspath_D_F(ch) - 1))^(-1);
                else
                    Dil_Eq(ch) = 1 + (Crosspath_I(ch) / Crosspath_d(ch))^(Crosspath_D_F(ch) - 1);
                end
            end
        else
            Not_Resolved_Count = Not_Resolved_Count + 1;
        end
    end

    MRCA_Locs_Reduced(~any(MRCA_Locs_Reduced, 2), :) = [];
    Rand_Combo_LoCs_Reduced = Rand_Combos_LoCs(Dil_Eq ~= 0, :);
    MRCA_Values_Reduced = MRCA_Values(Dil_Eq ~= 0);
    Crosspath_H_Reduced = Crosspath_H(Dil_Eq ~= 0);
    Crosspath_I_Reduced = Crosspath_I(Dil_Eq ~= 0);
    Crosspath_d_Reduced = Crosspath_d(Dil_Eq ~= 0);
    Crosspath_D_F = Crosspath_D_F(Crosspath_D_F ~= 0);
    Dil_Eq = Dil_Eq(Dil_Eq ~= 0);

    Dil_Eq_Real_Sols_Locs = find(imag(Dil_Eq) == 0);
    Dil_Eq_Real_Sols = Dil_Eq(Dil_Eq_Real_Sols_Locs);
    Dil_Eq_Complex_Sols_Locs = find(imag(Dil_Eq) ~= 0);
    Dil_Eq_Complex_Sols = Dil_Eq(Dil_Eq_Complex_Sols_Locs);

    Stat_Imag_Count(u) = numel(Dil_Eq_Complex_Sols);
    Stat_Real_Count(u) = numel(Dil_Eq_Real_Sols);
    Stat_Not_Resolved(u) = Not_Resolved_Count;
    Stat_Indeterminant(u) = Indeterminant_Count;
    Stat_Coherent_Count(u) = sum(Crosspath_d > 0);
    Stat_Divergent_Count(u) = sum(Crosspath_d < 0);
    Stat_Convergent_Count(u) = size(Rand_Combos_LoCs, 1) - Not_Resolved_Count;

    Stat_Percent_Coherent(u) = Stat_Coherent_Count(u) / (Stat_Coherent_Count(u) + Stat_Divergent_Count(u));
    Stat_Percent_Divergent(u) = Stat_Divergent_Count(u) / (Stat_Coherent_Count(u) + Stat_Divergent_Count(u));
    Stat_Percent_Real(u) = Stat_Real_Count(u) / (Stat_Real_Count(u) + Stat_Imag_Count(u));
    Stat_Percent_Imag(u) = Stat_Imag_Count(u) / (Stat_Real_Count(u) + Stat_Imag_Count(u));
    Stat_Percent_Convergent(u) = Stat_Convergent_Count(u) / size(Rand_Combos_LoCs, 1);
    Stat_Percent_Not_Resolved(u) = Stat_Not_Resolved(u) / size(Rand_Combos_LoCs, 1);
    Stat_Percent_Indeterminant(u) = Stat_Indeterminant(u) / size(Rand_Combos_LoCs, 1);

    % Save progress to a .mat file
    save('progress.mat', 'u', 'Stat_Percent_Coherent', 'Stat_Percent_Divergent', 'Stat_Percent_Real', ...
        'Stat_Percent_Imag', 'Stat_Percent_Convergent', 'Stat_Percent_Not_Resolved', 'Stat_Percent_Indeterminant', ...
        'Stat_Coherent_Count', 'Stat_Divergent_Count', 'Stat_Real_Count', 'Stat_Imag_Count', ...
        'Stat_Convergent_Count', 'Stat_Not_Resolved', 'Stat_Indeterminant');
end

% Calculate the average and standard deviation of the percentages
Stat_Avg_Percent_Coherent = mean(Stat_Percent_Coherent);
Stat_Std_Percent_Coherent = std(Stat_Percent_Coherent);
Stat_Avg_Percent_Divergent = mean(Stat_Percent_Divergent);
Stat_Std_Percent_Divergent = std(Stat_Percent_Divergent);
Stat_Avg_Percent_Real = mean(Stat_Percent_Real);
Stat_Std_Percent_Real = std(Stat_Percent_Real);
Stat_Avg_Percent_Imag = mean(Stat_Percent_Imag);
Stat_Std_Percent_Imag = std(Stat_Percent_Imag);
Stat_Avg_Percent_Convergent = mean(Stat_Percent_Convergent);
Stat_Std_Percent_Convergent = std(Stat_Percent_Convergent);
Stat_Avg_Percent_Not_Resolved = mean(Stat_Percent_Not_Resolved);
Stat_Std_Percent_Not_Resolved = std(Stat_Percent_Not_Resolved);
Stat_Avg_Percent_Indeterminant = mean(Stat_Percent_Indeterminant);
Stat_Std_Percent_Indeterminant = std(Stat_Percent_Indeterminant);