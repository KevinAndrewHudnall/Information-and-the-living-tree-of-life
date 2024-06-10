% BuildMultifractalTree

% This script builds a multifractal tree by implementing the iterated
% function system. The main objects returned are S, a matrix of scales
% representing all paths in the tree of life; and P, a matrix indicating
% the progeny had by each member at each iteration. Both P and S will have
% as many rows as there are iterates of the system. P will have as many
% columns as there are paths in the second to last iterate since the final
% iterate does not birth any children.

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%

% The only parameters are the maximum number of offspring and the maximum
% number of generations for the Galton-Watson branching process, as well as
% total number of system iterations. The structure blows up rapidly,
% setting practical limitations on the parameters.

MaxOffspring = 3; % Max number of offspring for every generated random tree
MaxGens = 2; % Max number of generations for every random tree
ITERATIONS = 15; % Total number of system iterations

%%%%%%%%%%%%%%%%%%%%% LEVEL 1: generate the first tree %%%%%%%%%%%%%%%%%%%%

Scale = rand(1, 1); % Randomly select a scale to generate tree
T = MakeRandomTree(MaxOffspring, MaxGens); % Generate random tree
Leaves = cell(ITERATIONS, 1); % Initialize cell array to hold leaf scales
Progenies = cell(ITERATIONS, 1); % Initialize cell array to hold progenies
Leaves{1} = ones(1, T) * Scale; % Store the scale of the leaves in tree 1
Progenies{1} = T; % Store the number of progeny in tree 1

%%%%%%%%%%%%%%% LEVEL 2 - ITERATIONS: generate all other trees %%%%%%%%%%%%

n = 1; % count iterations
while (n < ITERATIONS)
    for i = 1:size(Leaves{n}, 2) % loop over the leaves in the nth iterate
        
        % Randomly select a scale for each leaf in the nth iterate
        Scale = Leaves{n}(i) * rand(1, 1);
        
        % Generate a random tree for each one of these scales
        T = MakeRandomTree(MaxOffspring, MaxGens);
        
        % Store the leaves of the newly generated trees as iterate (n + 1)
        Leaves{n + 1} = [Leaves{n + 1}, ones(1, T) * Scale];
        
        % Update the number of progeny along each path as iterate (n + 1)
        Progenies{n + 1} = [Progenies{n + 1}, T];

    end
    n = n + 1;
end

%%%%%%%%%%% Organize Progenies and Leaves into matrices P and S %%%%%%%%%%%

% The number of offspring of every leaf at every iterate has been 
% recorded in Progenies, and likewise the scale of every leaf at every 
% iterate has been recorded in Leaves. These data structures need to be 
% converted to matrices that have N rows and as many columns as there are
% leaves in the final iterate. Doing so will make ancestral relationships
% explicit. Columns of the matrices will then correspond to distinct paths
% through the tree of life. Constructing these matrices requires
% determining how many times each member of every cell except the last of
% Progenies (and Leaves) must be repeated and populated in the matrix.
% Basically, the ancestral relationships are being back constructed using
% the information in Progenies, and then P and S are constructed via these
% ancestral relationships.

% tally holds the number of times each member of each cell of Progenies
% must be repeated in order to form the matrix P
tally = cell(size(Progenies, 1), 1);
for i = 1:size(Progenies, 1) % Loop over all rows in Progenies
    for j = 1:size(Progenies{i}, 2) % Loop over all columns in Progenies
        
        % Count where we start column-wise in our repititions
        count_1 = sum(Progenies{i}(1:(j - 1))) + 1; 
        
        % Count where we stop column-wise in our repititions relative to
        % our start location
        count_2 = sum(Progenies{i}(1:j));
        
        % Count the total number of repitions we need; initialize to 1
        count_3 = 1;
        
        clear value{i}
        
        % Loop over all cells of Progenies in iterate (i + 1)
        for k = (i + 1):(size(Progenies, 1))

            if (j == 1) % If we're on the first column
                
                % The total number of repeats equals the stop location
                count_3 = count_2;
                
                % Determine the next start-relative stop location
                count_2 = sum(Progenies{k}(count_1:count_2));
                
                % Start location is 1
                count_1 = 1;
                
            else
                % Number of repeats equals the start-relative stop location
                count_3 = count_2; 
                
                temp = count_1; % Set a temporary variable 
                
                if (k < size(Progenies, 1)) % If we're not at the last row
                    
                    % Update start location
                    count_1 = value{i}(1, k + 1) + 1; 
                    
                    % Update start-relative stop location
                    count_2 = count_1 + sum(Progenies{k}(temp:count_2))...
                        - 1;
                end
            end
            % Store the repeats required in a cell array
            value{i}(1, k) = count_3; 
        end
        % Store the repeats required in a cell array
        tally{i}(1, j) = count_3; 
    end
end

% Now use tally to build P
P = zeros(size(Progenies, 1), size(Progenies{end}, 2)); % Initialize P

% Generate the first row of P
P(1, 1:size(Progenies{end}, 2)) = ones(1, size(Progenies{end}, 2)) *...
    Progenies{1}(1, 1);

% Generate rows 2 through (last - 1) of P
for i = 2:(size(Progenies, 1) - 1)
    Continue = [];
    Begin = ones(1, tally{i}(1)) * Progenies{i}(1, 1);
    for j = 2:size(tally{i}, 2)
        
        Continue = [Continue, ones(1, (tally{i}(j) - tally{i}(j - 1))) *...
            Progenies{i}(1, j)];
        
    end
    P(i, 1:size(Progenies{end}, 2)) = [Begin, Continue]; 
end

% Generate the last row of P
P(end, 1:size(Progenies{end}, 2)) = Progenies{end}(1, :);

% Use Progenies to build S
S = cell(size(Leaves, 1), 1);
Temp_1 = Leaves;
for p = 0:(size(Leaves, 1) - 2)
    Temp_2 = cell(size(Leaves, 1) - p, 1);
    for i = (size(Leaves, 1) - p):-1:2
        for j = 1:size(Progenies{i + p, 1}, 2)

            Temp_2{i - 1} =...
                [Temp_2{i - 1}, ones(1, Progenies{i + p, 1}(j)) *...
                Temp_1{i - 1, 1}(j)];
        end
    end
    Temp_1 = Temp_2;
    
    S{(size(Leaves, 1) - (p + 1)), 1} =...
        Temp_2{(size(Leaves, 1) - (p + 1)), 1};
    
    clear Temp_2
end
S{size(Leaves, 1), 1} = Leaves{end, 1};
S = cell2mat(S);

% Every leaf in a final subtree is at the same scale, so lump them together
% to eliminate duplicates
S = unique(S', 'rows', 'stable')';
