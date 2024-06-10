function [Ancestor, Ancestor_Loc] = FindCommonAncestor(Leaf_1, Leaf_2)
% PRECONDITION:
%               "Leaf_1" and "Leaf_2" are any two members of the system
%               (i.e., entries in matrix S. They can be in the same iterate
%               or in different iterates.
%
% POSTCONDITION:
%               Their most recent common ancestor (MRC) "Ancestor" is
%               returned as well as its row and column location in S,
%               returned as "Ancestor_Loc". Zero is returned for "Ancestor"
%               if one of the leaves is the ancestor of all life since an
%               ancestor does not exist in this case.

global S % S is made global so it does not have to be passed to the 
         % function each time.

% Get the coordinates of the leaves in S
[row_coord_1, col_coord_1] = find(S == Leaf_1);
[row_coord_2, col_coord_2] = find(S == Leaf_2);

% Check if one of the leaves is the ancestor of all life. In this case a 
% parent does not exist
if ((row_coord_1(1) == 1) || (row_coord_2(1) == 1))
    Ancestor = 0;
    Ancestor_Loc(1) = 1;
    Ancestor_Loc(2) = 1;
    
else % Determine the parent of each leaf
    Parent_1 = S(row_coord_1(1) - 1, col_coord_1(1));
    Parent_2 = S(row_coord_2(1) - 1, col_coord_2(1));

    % Check if the parents are in the same subtree. If so, the MRC of the
    % leaves is their grandparent
    if (Parent_1 == Parent_2)

        Ancestor = S(row_coord_1(1) - 2, col_coord_1(1));
        Ancestor_Loc(1) = row_coord_1(1) - 2;
        Ancestor_Loc(2) = col_coord_1(1);

    % Check if the leaves are in the same iterate (i.e., the same row in S)
    elseif (row_coord_1(1) == row_coord_2(1))

        % Step backwards through the rows of S along the column of each
        % leaf to find their MRC as the first entry that is the same for
        % both columns.
        Ancestor = [];
        for i = (row_coord_1(1) - 1):-1:1
            if (S(i, col_coord_1(1)) == S(i, col_coord_2(1)))

                Ancestor = S(i, col_coord_1(1));
                Ancestor_Loc(1) = i;
                Ancestor_Loc(2) = col_coord_1(1);
                break
            end    
        end

    % If the leaves are not in the same iterate and Leaf 2 is in an earlier
    % iterate
    elseif (row_coord_1(1) > row_coord_2(1))

        % If the two leaves are in the same column (i.e., Leaf_2 is an
        % ancestor of Leaf_1
        if (ismember(col_coord_1(1), col_coord_2)) 

            Ancestor = Parent_2;
            Ancestor_Loc(1) = row_coord_2(1) - 1;
            Ancestor_Loc(2) = col_coord_2(1);

        % Else the two leaves are in different columns, so step backwards
        % through S to find their MRC as the first entry that is the same
        % for both columns.
        else
            Ancestor = [];
            for i = (row_coord_1(1)):-1:1
                if ((S(i, col_coord_1(1)) == S(i, col_coord_2(1))))

                    Ancestor = S(i, col_coord_2(1));
                    Ancestor_Loc(1) = i;
                    Ancestor_Loc(2) = col_coord_2(1);
                    break
                end    
            end
        end

    % If the leaves are not in the same iterate and Leaf 1 is in an earlier
    % iterate
    else %(row_coord_2(1) > row_coord_1(1))

        % If the two leaves are in the same column (i.e., Leaf_1 is an
        % ancestor of Leaf_2)
        if (ismember(col_coord_2(1), col_coord_1)) 

            Ancestor = Parent_1;
            Ancestor_Loc(1) = row_coord_1(1) - 1;
            Ancestor_Loc(2) = col_coord_1(1);
            
        % Else the two leaves are in different columns and different
        % iterates and Leaf 1 is earlier. Step backwards through the rows
        % of S along the column of each leaf to find their MRC as the first
        % entry that is the same for both columns.
        else
            Ancestor = [];
            for i = (row_coord_2(1)):-1:1
                if ((S(i, col_coord_2(1)) == S(i, col_coord_1(1))))

                    Ancestor = S(i, col_coord_2(1));
                    Ancestor_Loc(1) = i;
                    Ancestor_Loc(2) = col_coord_2(1);
                    break
                end    
            end
        end
    end
end
end