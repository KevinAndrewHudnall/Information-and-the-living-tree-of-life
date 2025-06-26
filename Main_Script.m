% Main_Script

% This script generates a large multifractal tree and computes all pairwise
% crosspath information quantities in parallel.
%
% This script performs the following steps:
%   1. Constructs a random multifractal tree using a random iterated function system (RIFS)
%   2. Calculates the fractal (correlation) dimension D_F for each lineage path
%   3. Determines where D_F has sufficiently converged to permit meaningful crosspath comparison
%   4. Chunks the data for scalable parallel processing using `parfeval`
%   5. Computes pairwise information-theoretic quantities (H, I, d), common ancestors,
%      and dilation equation solutions via `GetCrosspathQuantities`
%   6. Saves each chunk of results to disk to preserve memory and enable distributed analysis
%
% PARAMETERS:
%   ITERATIONS     - Number of recursive RIFS iterations (tree depth)
%   MaxOffspring   - Maximum number of offspring per node
%   MaxGens        - Number of generations per subtree
%   Chunk_Size     - Number of leaves per chunk for crosspath comparison
%   Tolerance      - Convergence threshold for accepting a comparison
%
% OUTPUT:
%   Each chunk is saved to disk as a `.mat` file in the specified directory.
%
% This script underlies the large-scale simulation pipeline for computing
% the full observer-relative structure of scale, entropy, and time used in
% Figures 6–10 of the manuscript *Information and the Living Tree of Life*.
%
% NOTE: This script is computationally intensive and memory-heavy.
% It may require multiple CPU cores and substantial RAM to run efficiently.
% Avoid running this script on laptops or limited-resource environments.
% Reducing ITERATIONS and/or MaxOffspring and MaxGens will reduce
% computation. It is advised to begin with a modest ITERATION size.

% Parameters used in manuscript.
ITERATIONS = 30;
MaxOffspring = 3;
MaxGens = 2;

% Build the random tree and get the scale matrix S, the offspring
% matrix P, and the entropy matrix H.
[S, P, H] = BuildMultifractalTreeFn(MaxOffspring, MaxGens, ITERATIONS);

% Get the fractal dimensions
[D_F, D_F_R_Squared, D_F_Avg_R2_Vector, D_F_Std_Devs, D_F_Random_Member, D_F_C_ep_Rand_Member,...
    D_F_ep_Rand_Member, D_F_Rand_Member, D_F_R_Squared_Rand_Member] = CalculateFractalDimsFn(S);

% Set a convergence tolerance
Tolerance = 0.1; % Fractal dimensions along a path must have converged to
                 % within 10% the final value of the path

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
        if(temp(i, j) == 0)
            temp(1:i, j) = 0;
        end
    end
end
Conv_Tol = temp;
clear temp Within_Tolerance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get all crosspath scale comparisons of leaves in the final iterate %
% and also record their locations in the original matrix S           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chunk the data so it can be parallelized; adjust based on computer specs
% and preference
Chunk_Size = 400;
Num_Chunks = floor(size(S, 2) / Chunk_Size) + 1;
Max_Concurrent_Chunk_Data = 1;     % maximize parallelization by setting it 
                                   % to the number of cpu cores

% Process chunks in batches
for Batch_Start = 1:Max_Concurrent_Chunk_Data:Num_Chunks
    Batch_End = min(Batch_Start + Max_Concurrent_Chunk_Data - 1, Num_Chunks);
    Num_Chunk_Data = Batch_End - Batch_Start + 1;

    % Initialize future objects for the current batch
    Chunk_Data(Num_Chunk_Data) = parallel.FevalFuture;
    for k = 1:Num_Chunk_Data
        Chunk_Data(k) = parfeval(@GetCrosspathQuantities, 1, S, D_F, ...
            Conv_Tol, Chunk_Size, Batch_Start + k - 1);
    end

    % Collect results as they complete
    for k = 1:Num_Chunk_Data
        [Completed_Index, Completed_Data] = fetchNext(Chunk_Data);

        Absolute_Chunk_Index = Batch_Start + Completed_Index - 1;
    
        % Save result to disk to free up RAM
        Dir_Path = 'H:\Data_Chunks\temp';
        % Use the absolute index relative to all chunks
        Filename = sprintf('Chunk_Struct_%d.mat', Absolute_Chunk_Index);
        Full_Path = fullfile(Dir_Path, Filename);
        save(Full_Path, 'Completed_Data');

        % Clear the completed future to save memory
        clear Completed_Data;
    end
end