% This script generates the full suite of figures from the manuscript
% *Information and the Living Tree of Life*, including plots of:
%   - Nested scale structure (Figure 5)
%   - Pathwise and crosspath entropy (Figure 7)
%   - Fractal dimensions (Figures 6, D1)
%   - Mutual information and joint entropy (Figure 8)
%   - Information distance (Figure 9)
%   - Dilation equation solutions (Figure 10)
%   - Observer-relative rates of time (Figure 12)
%   - Multifractal detrended fluctuation analysis (Figure E1)
%
% The script assumes that variables such as `Randomly_Reduced_S`,
% `D_F`, `Dil_Eq_Real_Sols`, `Generalized_Hurst_values`, etc., are
% preloaded into the workspace or available from earlier pipeline steps.
%
% Each figure is plotted using standard MATLAB graphics with log-scaling
% where appropriate. The script includes visualization of both raw and
% summary statistics across paths and pairwise comparisons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  SCALE SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Leaf_Scale_Matrix S, the paths Figure 5 TOP in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1) 
hold on 
for i = 1:size(Randomly_Reduced_S, 2)
    plot(1:size(Randomly_Reduced_S, 1), Randomly_Reduced_S(1:end, i), 'LineWidth', 1);
end
set(gca, 'yscale', 'log');
title('Scale random variable \beta');
xlabel('Iteration');
ylabel('scale \beta [form]');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Scale of leaves in final iterate as a line Figure 5 BOTTOM in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
plot(1:size(Randomly_Reduced_S, 2), Randomly_Reduced_S(end, :));
title('Scale of leaves in final iterate');
xlabel('Leaf number');
ylabel('Scale \beta [form]');
set(gca, 'yscale', 'log');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ENTROPY SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the path entropies Figure 7 TOP in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on 
for i = 1:size(Randomly_Reduced_H, 2)
    plot(1:size(Randomly_Reduced_H, 1), Randomly_Reduced_H(1:end, i));
end
title('Pathwise entropy of random variable \beta');
xlabel('Iteration');
ylabel('Path entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Entropies of leaves in final iterate as a line Figure 7 BOTTOM in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on
plot(1:size(Randomly_Reduced_H, 2), Randomly_Reduced_H(end, :));
title('Entropies of leaves in final iterate');
xlabel('Leaf number');
ylabel('Entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Fractal Dimensions Section   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the fractal dimension of every path at the final iterate Figure 6 TOP in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
plot(1:size(D_F, 2), D_F(end, :));
title('Fractal dimension of every path');
xlabel('Leaf number');
ylabel('D_F');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical plots for correlation dimension Figure D1 in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the regression for the randomly chosen path from script
% "CalculatFractalDims.m"
figure(6)
hold on
scatter(log(ep_Rand_Member), log(C_ep_Rand_Member), 20, 'filled');
% Get coefficients of a linear fit to the data.
Coeffs = polyfit(log(ep_Rand_Member), log(C_ep_Rand_Member), 1);
% Create x-axis
xFit = linspace(min(log(ep_Rand_Member)), max(log(ep_Rand_Member)), 1000);
% Get the estimated yFit value.
yFit = polyval(Coeffs , xFit);
plot(xFit, yFit, '--', 'LineWidth', 2); % Plot fitted line.
title('Regression to determine D_F for a randomly chosen path')
xlabel('log(\epsilon)');
ylabel('log(C(\epsilon))');

% Plot a histogram of the coefficients of determination at the final
% iterate
figure(7)
histogram(R_Squared(end, :));
xlabel('R^2')

% Plot the average coefficient of determination at each iterate
figure(8)
errorbar(3:size(Avg_R2_Vector, 2), Avg_R2_Vector(3:end),...
    Std_Devs(3:end), Std_Devs(3:end))  
title('Average coefficient of determination at each iterate');
xlabel('Number of points N');
ylabel('Average R^2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% CROSSPATH QUANTITIES SECTION  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot crosspath joint entropy combinations of leaves in final iterate Figure 8 LEFT in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
scatter(1:size(Crosspath_H_Reduced, 1), Crosspath_H_Reduced(:, 1), 1, 'filled')
set(gca, 'yscale', 'log');
title('Pairwise joint entropy of leaves in the final iterate');
xlabel('Pairwise Leaf comparison');
ylabel('Joint entropy [nats]');
axis tight
xlim([0 size(Crosspath_H_Reduced, 1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot crosspath mutual information combinations of leaves in final  %
% iterate Figure 8 RIGHT in manuscript                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
hold on
scatter(1:size(Crosspath_I_Reduced, 1), Crosspath_I_Reduced(:, 1), 1, 'filled');
set(gca, 'yscale', 'log');
title('Pairwise mutual information of leaves in final iterate');
xlabel('Pairwise leaf comparison');
ylabel('Mutual information [nats]');
axis tight
xlim([0 size(Crosspath_I_Reduced, 1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the information distance for combinations of leaves in the final   %
% iterate. Figure 9 in manuscript                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since some are positive and some are negative, and we'd like to plot both
% on a log scale, split Crosspath_d_Cells into positive and negative and
% keep spacing on the x-axis
for i = 1:size(Crosspath_d_Reduced, 1)
    if (Crosspath_d_Reduced(i, 1) > 0)

        Crosspath_d_Reduced_Pos(i, 1) = Crosspath_d_Reduced(i, 1);
    else
        Crosspath_d_Reduced_Pos(i, 1) = 0;
    end
end

for i = 1:size(Crosspath_d_Reduced, 1)
    if (Crosspath_d_Reduced(i, 1) < 0)

        Crosspath_d_Reduced_Neg(i, 1) = Crosspath_d_Reduced(i, 1);
    else
        Crosspath_d_Reduced_Neg(i, 1) = 0;
    end
end

% Plot the negative solutions on log scale
figure(11)
scatter(1:size(Crosspath_d_Reduced_Neg, 1), Crosspath_d_Reduced_Neg, 1, 'filled');
title('Negative crosspath information distance of leaves in final iterate');
xlabel('Pairwise leaf comparison');
ylabel('Information distance [nats]');
set(gca, 'yscale', 'log');
axis tight

% Plot the positive solutions on a log scale
figure(12)
scatter(1:size(Crosspath_d_Reduced_Pos, 1), Crosspath_d_Reduced_Pos, 1, 'filled');
title('Positive crosspath information distance of leaves in final iterate');
xlabel('Pairwise leaf comparison');
ylabel('Information distance [nats]');
set(gca, 'yscale', 'log');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the crosspath fractal dimensions Figure 6 BOTTOM in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)
scatter(1:size(Crosspath_D_F_Reduced, 1), Crosspath_D_F_Reduced, 5, 'filled');
title('Crosspath fractal dimension');
xlabel('Leaf to leaf comparison');
ylabel('Crosspath D_F');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Solutions to the scale dilation equation Figure 10 LEFT and RIGHT in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the real only solutions:
figure(14)
hold on
scatter(1:size(Dil_Eq_Real_Sols, 1),...
    Dil_Eq_Real_Sols, 3, 'filled');
scatter(1:size(Dil_Eq_Real_Sols_Inv, 1),...
    Dil_Eq_Real_Sols_Inv, 3, 'filled');
title('Real part of the dilation equation');
xlabel('Leaf comparison');
ylabel('Re(\Delta t^i_j)');
set(gca, 'yscale', 'log');

% Plot the complex without real only solutions:
figure(15)
hold on
scatter(real(Dil_Eq_Complex_Sols),...
    imag(Dil_Eq_Complex_Sols), 3, 'filled');
scatter(real(Dil_Eq_Complex_Sols_Inv),...
    imag(Dil_Eq_Complex_Sols_Inv), 3, 'filled');
title('Complex solutions to the dilation equation');
xlabel('Re(\Delta t^i_j)');
ylabel('Im(\Delta t^i_j)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTERED SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Rates of time elapse Figure 12 in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(16)
hold on
scatter(1:size(Rooted_Crosspath_Dil_Eq_Real, 1),...
    Rooted_Crosspath_Dil_Eq_Real, 5, 'filled');
title('Rates of time elapse');
xlabel('Pairwise comparison');
ylabel('Relative rate of time elapse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the DFA exponents Figure E1 in manuscript %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numRows = size(Generalized_Hurst_values, 1);  % Get the total number of rows in the matrix
% Generate 100,000 unique random row indices
sampleIndices = randperm(numRows, 50000);
% Use the indices to sample the rows from the matrix
sampledMatrix = Generalized_Hurst_values(sampleIndices, :);

figure(17)
hold on
for i = 1:size(sampledMatrix, 1)
    plot(1:size(sampledMatrix, 2), sampledMatrix(i, :));
end
title('MF-DFA');
xlabel('Moment q');
ylabel('Generalized Hurst values');
