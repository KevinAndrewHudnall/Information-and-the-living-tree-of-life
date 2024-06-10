% GetRootedQuantites

% This script randomly selects a leaf in the final iterate of S to be the
% root. It then gets the restriction of the dilation equation to only those
% solutions involving the root.

% Randomly select a root
Root_Location = randi(size(S, 2));
Root = S(end, Root_Location);
    
% Get the rooted dilation equation
count = 1;
for i = 1:size(Dilation_Eq, 1)
    
    if (S(ITERATIONS, Dilation_Eq_Locs_In_S(i, 1)) == Root ||...
            S(ITERATIONS, Dilation_Eq_Locs_In_S(i, 2)) == Root)
        
        Rooted_Dilation_Eq(count, 1) = Dilation_Eq(i);
        count = count + 1;
    end
end
clear count;

% Get the processed results
Rooted_Dilation_Eq_Inv = Rooted_Dilation_Eq.^(-1);

Rooted_Dilation_Eq_Locs = find(Dilation_Eq_Locs_In_S(:, 1) ==...
    Root_Location | Dilation_Eq_Locs_In_S(:, 2) == Root_Location);

Rooted_Locs_In_S = Dilation_Eq_Locs_In_S(Rooted_Dilation_Eq_Locs, :);

Rooted_Dilation_Eq_Complex_Sols =...
    Rooted_Dilation_Eq(imag(Rooted_Dilation_Eq) ~= 0);

Rooted_Dilation_Eq_Real_Sols =...
    Rooted_Dilation_Eq(imag(Rooted_Dilation_Eq) == 0);

Rooted_Dilation_Eq_Complex_Sols_Inv =...
    Rooted_Dilation_Eq_Inv(imag(Rooted_Dilation_Eq_Inv) ~= 0);

Rooted_Dilation_Eq_Real_Sols_Inv =...
    Rooted_Dilation_Eq_Inv(imag(Rooted_Dilation_Eq_Inv) == 0);