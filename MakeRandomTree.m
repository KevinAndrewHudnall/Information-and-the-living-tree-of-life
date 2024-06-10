% This function uses the Galton-Watson branching process to generate a
% random tree.

function [Tree] = MakeRandomTree(Max_Offspring, Max_Gens)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRECONDITION:                                                           %
%                Max_Offspring: Max number of offspring possible          %
%                Max_Gens: Max number of generations possible             %
% POSTCONDITION:                                                          %
%                Tree: The number of leaves generated                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Done = 0;
while(Done == 0)
    
    % Generate offspring probabilies
    Probs = rand(1, Max_Offspring + 1); 
    
    Cum_Probs = [cumsum(Probs) 1]; %Cumulative probabilities

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the branching %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    Parent = 0;
    Gens = randi(Max_Gens); % Determine number of generations
    Child_Count = 1;
    Next_Gen = cell(1, Gens);
    Gens_To_Extinction = 0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Galton-Watson Branching.                                                % 
% Modified from Ingemar Kaj Raimundas Gaigalas (2024)                     %
% https://www.mathworks.com/matlabcentral/fileexchange/2516-random-trees  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:Gens % Branching at each generation
        
        % Determine the linear indices of parents in previous generation.
        Parent_Index = length(Parent) - Child_Count + 1:length(Parent);
        
        % Update Parent_Count to continue getting linear indices of parents
        Parent_Count = Child_Count;
        
        % Generate the offspring distribution
        Child_Dist = rand(1, Parent_Count);

        Child_Count = 0;
        for j = 1:Max_Offspring
            
            % Set Index to contain all values of Parent_Index that are
            % greater than the jth value of Cum_Probs but less than the
            % (j+1)th value of Cum_Probs.
            Index = Parent_Index((Child_Dist > Cum_Probs(j)) &...
                (Child_Dist <= Cum_Probs(j + 1)));

            if (~isempty(Index)) % If at least one parent birthed j kids
                
                % Get the number of offspring for each parent by making a
                % vector of Parent concatenated with the vector Index
                % repeated j times.
                Parent = [Parent repmat(Index, 1, j)];
                
                % Update Child_Count with additional children
                Child_Count = Child_Count + length(Index) * j;
            end
        end
        Next_Gen{i} = Parent;

        % If the tree has gone extinct, record the generations to
        % extinction and break from the loop.
        if (Child_Count == 0)
            Gens_To_Extinction = i;
            break;
        end
    end

    % If the tree didn't go extinct in the first generation
    if (Gens_To_Extinction ~= 1) 
                        
        %Determine leaves and set it to Tree
        Parent = rem(Parent + length(Parent), length(Parent) + 1) + 1;
        Isaleaf = ones(1, length(Parent) + 1);
        Isaleaf(Parent) = zeros(length(Parent), 1);
        
        [~, Tree] = size(Isaleaf(Isaleaf == 1));
        
        Done = 1;
    end
end

