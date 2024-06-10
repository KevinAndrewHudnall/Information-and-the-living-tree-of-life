%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Statistical_Confirmation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script confirms that the percent divergence versus percent
% convergence reported in the manuscript is statistically representative

global S;

% Decide how many systems to simulate to compare the statistics. If
% ITERATIONS is large, and/or MaxOffspring and MaxGens are large, then the
% systems generated will be large, presenting potential computational
% difficulties.
Number_Of_Cases = 5; 

for u = 1:Number_Of_Cases
    
    % Build the random tree
    BuildMultifractalTree;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Entropy at each level for each path %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H = zeros(size(S));
    for i = 1:size(S, 2)
        for j = 2:size(S, 1)
            H(j, i) = log(S(j - 1, i));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the information quantities for all pairwise comparisons %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Locs_In_S] = (unique(S(end, :), 'stable')); 
    Crosspath_Combos = nchoosek(S(end, :), 2);
    Locations_In_S = nchoosek(Locs_In_S(:, 1), 2);
    for i = 1:size(Crosspath_Combos, 1)
        Parent_1 = S(ITERATIONS - 1, Locations_In_S(i, 1));
        Parent_2 = S(ITERATIONS - 1, Locations_In_S(i, 2));
        [Common_Ancestor, Common_Ancestor_Loc] = ...
            FindCommonAncestor(Crosspath_Combos(i, 1), Crosspath_Combos(i, 2));
        Crosspath_H = ((2 * Parent_1 * Parent_2) / Common_Ancestor^2) *...
            log(Common_Ancestor);
        Crosspath_I = (Parent_1 * Parent_2 / Common_Ancestor^2) *...
            (log((Parent_1 * Parent_2) / Common_Ancestor^2));
        Crosspath_d(i, 1) = Crosspath_H - Crosspath_I;
    end
    
    % Save the results for statistical analysis
    Crosspath_d_For_Stat{u} = Crosspath_d;
    Coherent_For_Stat(u) = size(find(Crosspath_d > 0), 1) / size(Crosspath_d, 1);
    Divergent_For_Stat(u) = size(find(Crosspath_d < 0), 1) / size(Crosspath_d, 1);
end
clearvars -except Crosspath_d_For_Stat Coherent_For_Stat Divergent_For_Stat
