%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% CROSSPATH QUANTITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script calculates the crosspath information quantities for all 
% pairwise comparisons of leaves in the final iterate: joint entropy,
% mutual information, and information distance. It also calculates the
% dilation equations for all valid pairwise comparisons. A comparison is
% valid if it is within a user specified convergence tolerance. We used a
% tolerance of 50% for the results in the manuscript.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get all crosspath scale comparisons of leaves in the final iterate %
% and also record their locations in the original matrix S           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global S % Make S global to avoid passing it between functions

% Set up all pairwise comparisons
Locs_In_S = nchoosek(1:size(S, 2), 2);
S_Final = S(ITERATIONS, :);
Crosspath_Combos = S_Final(Locs_In_S);
clear S_Final;

Tolerance = 0.5; % Fractal dimensions along a path must have converged to
                 % within 50% the final value of the path

% Obtain a matrix that gives D_F for those within tolerance                 
Conv_Tol = zeros(size(D_F));
for i = 1:size(D_F, 2)
    for j  = 1:size(D_F, 1)
        Conv_Tol(j, i) = abs(D_F(end, i) - D_F(j, i));
    end
end

for i = 1:size(Conv_Tol, 1)
    for j = 1:size(Conv_Tol, 2)
        if(Conv_Tol(i, j) < (Tolerance * D_F(end, j)))
            temp(i, j) = D_F(i, j);
        else
            temp(i, j) = 0;
        end
    end
end

for i = size(temp, 1):-1:1
    for j = 1:size(temp, 2)
        if(temp(i, j) == 0)
            temp(1:i, j) = 0;
        end
    end
end
Conv_Tol = temp;
clear temp;
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the crosspath quantities for the pairwise comparison of %
% leaves in the final iterate                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize solution arrays
Crosspath_H = zeros(size(Crosspath_Combos, 1), 1);
Crosspath_I = zeros(size(Crosspath_Combos, 1), 1);
Crosspath_d = zeros(size(Crosspath_Combos, 1), 1);
Crosspath_D_F = zeros(size(Crosspath_Combos, 1), 1);
Dilation_Eq = zeros(size(Crosspath_Combos, 1), 1);
Dilation_Eq_Locs_In_S = zeros(size(Crosspath_Combos, 1), 1);
Not_Resolved = 0;
Indeterminant = 0;

for i = 1:size(Crosspath_Combos, 1)

    % Determine their parents
    Parent_1 = S(ITERATIONS, Locs_In_S(i, 1));
    Parent_2 = S(ITERATIONS, Locs_In_S(i, 2));

    % Determine their common ancestor and its location in S
    [Common_Ancestor, Common_Ancestor_Loc] = ...
        FindCommonAncestor(Crosspath_Combos(i, 1), Crosspath_Combos(i, 2));

    % Get the information quantities
    Crosspath_H(i, 1) = ((2 * Parent_1 * Parent_2) /...
        Common_Ancestor^2) * log(Common_Ancestor);
    Crosspath_I(i, 1) = (Parent_1 * Parent_2 / Common_Ancestor^2) *...
        (log((Parent_1 * Parent_2) / Common_Ancestor^2));
    Crosspath_d(i, 1) = Crosspath_H(i, 1) - Crosspath_I(i, 1);

    % If the most recent common ancestor is the ancestor of all life,
    % or the child of the ancestor of all life, these comparisons must be
    % eliminated since a fractal dimension does not exist in these cases.
    if (Common_Ancestor == S(1, 1) || ismember(Common_Ancestor, S(2, :)))

        Not_Resolved = Not_Resolved + 1;

    % If D_F for the most recent common ancestor has not converged to
    % within tolerance along the path for both entities being compared,
    % then these comparisons must be eliminated.    
    elseif (Conv_Tol(Common_Ancestor_Loc(1), Locs_In_S(i, 1)) == 0 ||...
            Conv_Tol(Common_Ancestor_Loc(1), Locs_In_S(i, 2)) == 0)

        Not_Resolved = Not_Resolved + 1;

    % If H[i] = H[j] no decision criterion exists for the dilation
    % equation so the comparison must be omitted.    
    elseif (abs(log(S(p - 1, Locs_In_S(i, 1))))) ==...
            abs(log(S(p - 1, Locs_In_S(i, 2))))

        Indeterminant = Indeterminant + 1;
        
    else % The comparison can be made

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get the joint fractal dimension  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Crosspath_D_F(i) = D_F(Common_Ancestor_Loc(1),...
            Common_Ancestor_Loc(2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the dilation equation for this pairwise comparisons. %                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % See which one has greater entropy and select either the dilation
        % equation or its inverse as the solution.
        if (abs(log(S(p - 1, Locs_In_S(i, 1))))) > ...
                abs(log(S(p - 1, Locs_In_S(i, 2))))

            Dilation_Eq(i) = (1 + (Crosspath_I(i) / Crosspath_d(i))^...
                (Crosspath_D_F(i) - 1))^(-1);
        else
            Dilation_Eq(i) =  1 + (Crosspath_I(i) / Crosspath_d(i))^...
                (Crosspath_D_F(i) - 1);
        end
        
        % Record their locations in S
        Dilation_Eq_Locs_In_S(i, 1) = Locs_In_S(i, 1);
        Dilation_Eq_Locs_In_S(i, 2) = Locs_In_S(i, 2);
    end
end

% Eliminate Indeterminant, Not_Resolved
Dilation_Eq = Dilation_Eq(Dilation_Eq ~= 0);
Crosspath_D_F = Crosspath_D_F(Crosspath_D_F ~= 0);

count = 1;
for i = 1:size(Dilation_Eq_Locs_In_S, 1)
    if (Dilation_Eq_Locs_In_S(i, 1) ~= 0)
        temp(count, :) = Dilation_Eq_Locs_In_S(i, :);
        count = count + 1;
    end
end
Dilation_Eq_Locs_In_S = temp;
clear temp

%%%%%%%%%%%%%%%%%%%%
% Organize results %
%%%%%%%%%%%%%%%%%%%%

% Determine percent in the positive regime = percent coherent evolution
Percent_Coherent = size(find(Crosspath_d > 0), 1) / size(Crosspath_d, 1);
% Determine percent in the negative regime = percent divergent evolution
Percent_Divergent = size(find(Crosspath_d < 0), 1) / size(Crosspath_d, 1);
% Get the REAL ONLY solutions 
Dilation_Eq_Real_Sols_Locs = find(imag(Dilation_Eq) == 0);
Dilation_Eq_Real_Sols = Dilation_Eq(Dilation_Eq_Real_Sols_Locs);
% Get the IMAG ONLY solutions
Dilation_Eq_Complex_Sols_Locs = find(imag(Dilation_Eq) ~= 0);
Dilation_Eq_Complex_Sols = Dilation_Eq(Dilation_Eq_Complex_Sols_Locs);
% Get the inverse solutions
Dilation_Eq_Inv = Dilation_Eq.^(-1);
Dilation_Eq_Real_Sols_Inv = Dilation_Eq_Real_Sols.^(-1);
Dilation_Eq_Complex_Sols_Inv = Dilation_Eq_Complex_Sols.^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the percent that are imaginary only and real only %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Percent_Imag = size(Dilation_Eq_Complex_Sols, 1) / size(Dilation_Eq, 1);
Percent_Real = size(Dilation_Eq_Real_Sols, 1) / size(Dilation_Eq, 1);