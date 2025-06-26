function [Data_Chunk] = GetCrosspathQuantities(S, D_F, Conv_Tol, Chunk_Size, k)

% This function computes all pairwise crosspath information
% quantities, common ancestors, and solutions to the dilation equation
% for a chunk of leaf-to-leaf comparisons in a multifractal tree.
%
% INPUTS:
%   S          - scale matrix from BuildMultifractalTreeFn
%   D_F        - fractal dimension matrix from CalculateFractalDimsFn
%   Conv_Tol   - convergence tolerance matrix to screen comparisons
%   Chunk_Size - number of leaves to include in each chunk
%   k          - index of the current chunk (1-based)
%
% OUTPUT:
%   Data_Chunk - struct containing:
%       • Pairwise comparison metadata (locations, parents, MRCAs)
%       • Crosspath entropy H, mutual information I, and information distance d
%       • Dilation equation solutions and their inverses (real/complex split)
%       • Crosspath fractal dimensions at MRCAs
%       • Counts of coherent (positive d), divergent (negative d),
%         imaginary, and indeterminant results
%
% The function iterates over all pairs of leaves within the current chunk,
% plus all prior chunks, to compute:
%   - MRCA scale values
%   - Entropy-based quantities (H, I, d)
%   - Dilation equation solutions based on entropy and D_F
%
% It filters comparisons using Conv_Tol and labels unresolved,
% indeterminant, and imaginary outcomes. This function is designed
% to be run in parallel across chunks for large-scale pairwise analysis.

Start_Idx = (k - 1) * Chunk_Size + 1;
End_Idx = k * Chunk_Size;

if (End_Idx > size(S, 2))
    End_Idx = size(S, 2);
end

Chunk_Combo_Locs = nchoosek(Start_Idx:End_Idx, 2);
% Chunk_Combo_Values = S_Final(Chunk_Combo_Locs);

for i = 1:(k - 1)
    % Define the previous chunk
    Prev_Start_Idx = (i - 1) * Chunk_Size + 1;

    Prev_End_Idx = min(i * Chunk_Size, length(S(end, :)));

    Prev_Chunk = Prev_Start_Idx:Prev_End_Idx;

    % Generate pairwise comparisons between current chunk and previous chunk
    Cross_Chunk_Pairs = [];
    for p = 1:length(Prev_Chunk)
        Cross_Chunk_Pairs = [Cross_Chunk_Pairs; [repmat(Prev_Chunk(p), ...
            length(Start_Idx:End_Idx), 1), (Start_Idx:End_Idx)']];
    end

    % Store the Cross_Chunk_Pairs
    Chunk_Combo_Locs = [Chunk_Combo_Locs; Cross_Chunk_Pairs];
end

% Get a Chunk_Size by 2 array of parents
Parent_1_Chunk = transpose(S(size(S, 1) - 1, Chunk_Combo_Locs(:, 1)));
Parent_2_Chunk = transpose(S(size(S, 1) - 1, Chunk_Combo_Locs(:, 2)));

% Get a Chunk_Size by 2 array of most recent common ancestors (MRCA)
MRCA_Chunk_Locs = zeros(size(Chunk_Combo_Locs, 1), 2);
MRCA_Chunk_Values = zeros(size(Chunk_Combo_Locs, 1), 1);
MRCA_Chunk_Locs_Reduced = zeros(size(Chunk_Combo_Locs, 1), 2);
Crosspath_H_Chunk = zeros(size(Chunk_Combo_Locs, 1), 1);
Crosspath_I_Chunk = zeros(size(Chunk_Combo_Locs, 1), 1);
Crosspath_d_Chunk = zeros(size(Chunk_Combo_Locs, 1), 1);
Crosspath_D_F_Chunk = zeros(size(Chunk_Combo_Locs, 1), 1);
Dil_Eq_Chunk = zeros(size(Chunk_Combo_Locs, 1), 1);
% Dil_Eq_Chunk_Locs = zeros(size(Chunk_Combo_Locs, 1), 2);
Not_Resolved = 0;
Indeterminant = 0;
for ch = 1:size(Chunk_Combo_Locs, 1)

    % Check if the parents are in the same subtree. If so, the MRC of the
    % leaves is their grandparent
    if (Parent_1_Chunk(ch, 1) == Parent_2_Chunk(ch, 1))

        Ancestors_Loc_New(1, 1) = size(S, 1) - 2;
        Ancestors_Loc_New(1, 2) = Chunk_Combo_Locs(ch, 1);
    else
        % Step backwards through the rows of S along the column of each
        % leaf to find their MRC as the first entry that is the same for
        % both columns.
        for i = (size(S, 1) - 1):-1:1
            if (S(i, Chunk_Combo_Locs(ch, 1)) == S(i, Chunk_Combo_Locs(ch, 2)))

                Ancestors_Loc_New(1, 1) = i;
                Ancestors_Loc_New(1, 2) = Chunk_Combo_Locs(ch, 1);
                break
            end    
        end
    end

    % Organize arrays of most recent common ancestors' locations and values
    MRCA_Chunk_Locs(ch, :) = Ancestors_Loc_New;
    MRCA_Chunk_Values(ch, 1) = S(MRCA_Chunk_Locs(ch, 1), MRCA_Chunk_Locs(ch, 2));

    % Calculate the information quantities
    Crosspath_H_Chunk(ch, 1) = ((2 * Parent_1_Chunk(ch, 1) * ...
        Parent_2_Chunk(ch, 1)) / MRCA_Chunk_Values(ch, 1)^2) *...
        log(MRCA_Chunk_Values(ch, 1));
    Crosspath_I_Chunk(ch, 1) = (Parent_1_Chunk(ch, 1) * Parent_2_Chunk(ch, 1) / ...
        MRCA_Chunk_Values(ch, 1)^2) * (log((Parent_1_Chunk(ch, 1) * ...
        Parent_2_Chunk(ch, 1)) / MRCA_Chunk_Values(ch, 1)^2));
    Crosspath_d_Chunk(ch, 1) = Crosspath_H_Chunk(ch, 1) - Crosspath_I_Chunk(ch, 1);

    % Get reduced arrays
    if (Conv_Tol(MRCA_Chunk_Locs(ch, 1), MRCA_Chunk_Locs(ch, 2)) ~= 0)

        MRCA_Chunk_Locs_Reduced(ch, :) = [MRCA_Chunk_Locs(ch, 1), MRCA_Chunk_Locs(ch, 2)];

        if (MRCA_Chunk_Locs_Reduced(ch, 1) == 1 || MRCA_Chunk_Locs_Reduced(ch, 1) == 2)

            Not_Resolved = Not_Resolved + 1;

        elseif (abs(log(S(size(S, 1) - 1, Chunk_Combo_Locs(ch, 1))))) ==...
                abs(log(S(size(S, 1) - 1, Chunk_Combo_Locs(ch, 2))))

            Indeterminant = Indeterminant + 1;

        elseif (MRCA_Chunk_Locs_Reduced(ch, 1) ~= 0 ...
                && MRCA_Chunk_Locs_Reduced(ch, 2) ~= 0) % The comparison can be made

            Crosspath_D_F_Chunk(ch) = D_F(MRCA_Chunk_Locs_Reduced(ch, 1),...
                MRCA_Chunk_Locs_Reduced(ch, 2));

            % See which one has greater entropy and select either the dilation
            % equation or its inverse as the solution.
            if (abs(log(S(size(S, 1) - 1, Chunk_Combo_Locs(ch, 1)))) > ...
                    abs(log(S(size(S, 1) - 1, Chunk_Combo_Locs(ch, 2)))))

                Dil_Eq_Chunk(ch) = (1 + (Crosspath_I_Chunk(ch, 1) / Crosspath_d_Chunk(ch, 1))^...
                    (Crosspath_D_F_Chunk(ch) - 1))^(-1);
            else
                Dil_Eq_Chunk(ch) =  1 + (Crosspath_I_Chunk(ch, 1) / Crosspath_d_Chunk(ch, 1))^...
                    (Crosspath_D_F_Chunk(ch) - 1);
            end
            % % Record their locations in S
            % Dil_Eq_Chunk_Locs(ch, 1) = Chunk_Combo_Locs(ch, 1);
            % Dil_Eq_Chunk_Locs(ch, 2) = Chunk_Combo_Locs(ch, 2);
        end
    else
        Not_Resolved = Not_Resolved + 1;
    end
end

% Reduce the number of pairwise comparisons to only those within
% convergence tolerance.
MRCA_Chunk_Locs_Reduced(~any(MRCA_Chunk_Locs_Reduced, 2), :) = [];
Chunk_Combo_Locs_Reduced = Chunk_Combo_Locs(Dil_Eq_Chunk(:, 1) ~= 0, :);
MRCA_Chunk_Values_Reduced = MRCA_Chunk_Values(Dil_Eq_Chunk ~= 0);
Crosspath_H_Chunk_Reduced = Crosspath_H_Chunk(Dil_Eq_Chunk ~= 0);
Crosspath_I_Chunk_Reduced = Crosspath_I_Chunk(Dil_Eq_Chunk ~= 0);
Crosspath_d_Chunk_Reduced = Crosspath_d_Chunk(Dil_Eq_Chunk ~= 0);
% Parent_1_Chunk_Reduced = Parent_1_Chunk(Parent_1_Chunk ~= 0);
% Parent_2_Chunk_Reduced = Parent_2_Chunk(Parent_2_Chunk ~= 0);
Crosspath_D_F_Chunk = Crosspath_D_F_Chunk(Crosspath_D_F_Chunk ~= 0);
Dil_Eq_Chunk = Dil_Eq_Chunk(Dil_Eq_Chunk ~= 0);
% Dil_Eq_Chunk_Locs(~any(Dil_Eq_Chunk_Locs, 2), :) = [];

%%%%%%%%%%%%%%%%%%%%
% Organize results %
%%%%%%%%%%%%%%%%%%%%

% Determine percent in the positive regime = percent coherent evolution
Coherent_Count = sum(Crosspath_d_Chunk > 0);
% Determine percent in the negative regime = percent divergent evolution
Divergent_Count = sum(Crosspath_d_Chunk < 0);

% Get the REAL ONLY solutions 
Dil_Eq_Chunk_Real_Sols_Locs = find(imag(Dil_Eq_Chunk) == 0);
Dil_Eq_Chunk_Real_Sols = Dil_Eq_Chunk(Dil_Eq_Chunk_Real_Sols_Locs);
% Get the IMAG ONLY solutions
Dil_Eq_Chunk_Complex_Sols_Locs = find(imag(Dil_Eq_Chunk) ~= 0);
Dil_Eq_Chunk_Complex_Sols = Dil_Eq_Chunk(Dil_Eq_Chunk_Complex_Sols_Locs);
% Get the inverse solutions
Dil_Eq_Chunk_Inv = Dil_Eq_Chunk.^(-1);
Dil_Eq_Chunk_Real_Sols_Inv = Dil_Eq_Chunk_Real_Sols.^(-1);
Dil_Eq_Chunk_Complex_Sols_Inv = Dil_Eq_Chunk_Complex_Sols.^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the percent that are imaginary only and real only %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imag_Count = size(Dil_Eq_Chunk_Complex_Sols, 1);
Real_Count = size(Dil_Eq_Chunk_Real_Sols, 1);

% Create the struct
Data_Chunk = struct();

Data_Chunk.Chunk_Combo_Locs = Chunk_Combo_Locs;
Data_Chunk.Parent_1_Chunk = Parent_1_Chunk;
Data_Chunk.Parent_2_Chunk = Parent_2_Chunk;
Data_Chunk.MRCA_Chunk_Locs = MRCA_Chunk_Locs;
Data_Chunk.MRCA_Chunk_Values = MRCA_Chunk_Values;
Data_Chunk.Crosspath_H_Chunk = Crosspath_H_Chunk;
Data_Chunk.Crosspath_I_Chunk = Crosspath_I_Chunk;
Data_Chunk.Crosspath_d_Chunk = Crosspath_d_Chunk;
% Data_Chunk.Chunk_Combo_Values = Chunk_Combo_Values;

Data_Chunk.Chunk_Combo_Locs_Reduced = Chunk_Combo_Locs_Reduced;
% Data_Chunk.Chunk_Combo_Values_Reduced
Data_Chunk.MRCA_Chunk_Locs_Reduced = MRCA_Chunk_Locs_Reduced;
Data_Chunk.MRCA_Chunk_Values_Reduced = MRCA_Chunk_Values_Reduced;
Data_Chunk.Crosspath_H_Chunk_Reduced = Crosspath_H_Chunk_Reduced;
Data_Chunk.Crosspath_I_Chunk_Reduced = Crosspath_I_Chunk_Reduced;
Data_Chunk.Crosspath_d_Chunk_Reduced = Crosspath_d_Chunk_Reduced;
% Data_Chunk.Parent_1_Chunk_Reduced = Parent_1_Chunk_Reduced;
% Data_Chunk.Parent_2_Chunk_Reduced = Parent_2_Chunk_Reduced;

Data_Chunk.Crosspath_D_F_Chunk = Crosspath_D_F_Chunk;
Data_Chunk.Dil_Eq_Chunk = Dil_Eq_Chunk;
% Data_Chunk.Dil_Eq_Chunk_Locs = Dil_Eq_Chunk_Locs;

Data_Chunk.Coherent_Count = Coherent_Count;
Data_Chunk.Divergent_Count = Divergent_Count;

Data_Chunk.Dil_Eq_Chunk_Real_Sols_Locs = Dil_Eq_Chunk_Real_Sols_Locs;
Data_Chunk.Dil_Eq_Chunk_Real_Sols = Dil_Eq_Chunk_Real_Sols;
Data_Chunk.Dil_Eq_Chunk_Complex_Sols_Locs = Dil_Eq_Chunk_Complex_Sols_Locs;
Data_Chunk.Dil_Eq_Chunk_Complex_Sols = Dil_Eq_Chunk_Complex_Sols;
Data_Chunk.Dil_Eq_Chunk_Inv = Dil_Eq_Chunk_Inv;
Data_Chunk.Dil_Eq_Chunk_Real_Sols_Inv = Dil_Eq_Chunk_Real_Sols_Inv;
Data_Chunk.Dil_Eq_Chunk_Complex_Sols_Inv = Dil_Eq_Chunk_Complex_Sols_Inv;
Data_Chunk.Imag_Count = Imag_Count;
Data_Chunk.Real_Count = Real_Count;

Data_Chunk.Not_Resolved = Not_Resolved;
if exist("Indeterminant", "var")
    Data_Chunk.Indeterminant = Indeterminant;
end