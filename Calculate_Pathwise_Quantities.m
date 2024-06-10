% This script calulates the relevant pathwise quantities: path entropies,
% systemwide entropy and per capita entropy, 

global S; % Make the scale matrix global

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%       Entropy H Section       %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Entropy at each level for each path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get H to make calculated fractal dimensions easier
H = zeros(size(S));
for i = 1:size(S, 2)
    for j = 2:size(S, 1)
        H(j, i) = log(S(j - 1, i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system-wide entropy values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine entropy of entire system
System_Wide_H_Vector = sum(H, 2);

% Determine the total number of leaves at each iterate
Total_Leaves_By_Iterate = zeros(ITERATIONS, 1);
for i = 1:ITERATIONS
    Total_Leaves_By_Iterate(i, 1) = size(Progenies{i, 1}, 2);
end

% Determine total system entropy averaged by leaves per iterate (i.e.,
% per capita)
Per_Capita_H_Vector = System_Wide_H_Vector(:) ./ Total_Leaves_By_Iterate;