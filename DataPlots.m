%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  SCALE SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Leaf_Scale_Matrix S, the paths %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on 
for i = 1:size(S, 2)
    plot(1:size(S, 1), S(1:end, i), 'LineWidth', 1);
end
set(gca, 'yscale', 'log');
title('Scale random variable \beta');
xlabel('Iteration');
ylabel('scale \beta [form]');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Scale of leaves in final iterate as a line %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
plot(1:size(S, 2), S(end, :));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the path entropies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on 
for i = 1:size(H, 2)
    plot(1:size(H, 1), H(1:end, i));
end
title('Pathwise entropy of random variable \beta');
xlabel('Iteration');
ylabel('Path entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Entropies of leaves in final iterate as a line %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on
plot(1:size(H, 2), H(end, :));
title('Entropies of leaves in final iterate');
xlabel('Leaf number');
ylabel('Entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot system wide H by iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
plot(1:size(System_Wide_H_Vector, 1), System_Wide_H_Vector(1:end))
title('System-wide entropy by iterate');
xlabel('Iterate');
ylabel('Entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot system wide per capita H by iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
plot(1:size(Per_Capita_H_Vector, 1), Per_Capita_H_Vector(1:end))
title('System-wide per capita entropy by iterate');
xlabel('Iterate');
ylabel('Entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Fractal Dimensions Section   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the fractal dimension of every path at the final iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
plot(1:size(D_F, 2), D_F(end, :));
title('Fractal dimension of every path');
xlabel('Leaf number');
ylabel('D_F');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical plots for correlation dimension %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the regression for the randomly chosen path from script
% "CalculatFractalDims.m"
figure(8)
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
figure(9)
histogram(R_Squared(end, :));
xlabel('R^2')

% Plot the average coefficient of determination for at each iterate
figure(10)
errorbar(1:size(Avg_R2_Vector(3:end), 2), Avg_R2_Vector(3:end),...
    Std_Devs(3:end), Std_Devs(3:end))  
title('Average coefficient of determination at each iterate');
xlabel('Number of points N');
ylabel('Average R^2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% CROSSPATH QUANTITIES SECTION  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                               %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot crosspath joint entropy combinations of leaves in final iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
scatter(1:size(Crosspath_H, 1), Crosspath_H(:, 1), 1, 'filled')
set(gca, 'yscale', 'log');
title('Pairwise joint entropy of leaves in the final iterate');
xlabel('Pairwise Leaf comparison');
ylabel('Joint entropy [nats]');
axis tight
xlim([0 size(Crosspath_H, 1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot crosspath mutual information combinations of leaves in final  %
% iterate                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12)
hold on
scatter(1:size(Crosspath_I, 1), Crosspath_I(:, 1), 1, 'filled');
set(gca, 'yscale', 'log');
title('Pairwise mutual information of leaves in final iterate');
xlabel('Pairwise leaf comparison');
ylabel('Mutual information [nats]');
axis tight
xlim([0 size(Crosspath_I, 1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot crosspath variation of information combinations of leaves in final %
% iterate                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since some are positive and some are negative, and we'd like to plot both
% on a log scale, split Crosspath_d_Cells into positive and negative and
% keep spacing on the x-axis
for i = 1:size(Crosspath_d, 1)
    if (Crosspath_d(i, 1) > 0)

        Crosspath_d_Pos(i, 1) = Crosspath_d(i, 1);
    else
        Crosspath_d_Pos(i, 1) = 0;
    end
end

for i = 1:size(Crosspath_d, 1)
    if (Crosspath_d(i, 1) < 0)

        Crosspath_d_Neg(i, 1) = Crosspath_d(i, 1);
    else
        Crosspath_d_Neg(i, 1) = 0;
    end
end

% Plot the negative solutions on log scale
figure(13)
scatter(1:size(Crosspath_d_Neg, 1), Crosspath_d_Neg, 1, 'filled');
title('Negative crosspath information distance of leaves in final iterate');
xlabel('Pairwise leaf comparison');
ylabel('Information distance [nats]');
set(gca, 'yscale', 'log');
axis tight

% Plot the positive solutions on a log scale
figure(14)
scatter(1:size(Crosspath_d_Pos, 1), Crosspath_d_Pos, 1, 'filled');
title('Positive crosspath information distance of leaves in final iterate');
xlabel('Pairwise leaf comparison');
ylabel('Information distance [nats]');
set(gca, 'yscale', 'log');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the crosspath fractal dimensions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(15)
scatter(1:size(Crosspath_D_F, 1), Crosspath_D_F, 5, 'filled');
title('Crosspath fractal dimension');
xlabel('Leaf to leaf comparison');
ylabel('Crosspath D_F');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Solutions to the scale dilation equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the real only solutions:
figure(16)
hold on
scatter(1:size(Dilation_Eq_Real_Sols, 1),...
    Dilation_Eq_Real_Sols, 3, 'filled');
scatter(1:size(Dilation_Eq_Real_Sols_Inv, 1),...
    Dilation_Eq_Real_Sols_Inv, 3, 'filled');
title('Real part of the dilation equation');
xlabel('Leaf comparison');
ylabel('Re(\Delta t^i_j)');
set(gca, 'yscale', 'log');

% Plot the complex without real only solutions:
figure(17)
hold on
scatter(real(Dilation_Eq_Complex_Sols),...
    imag(Dilation_Eq_Complex_Sols), 3, 'filled');
scatter(real(Dilation_Eq_Complex_Sols),...
    imag(Dilation_Eq_Complex_Sols), 3, 'filled');
title('Complex solutions to the dilation equation');
xlabel('Re(\Delta t^i_j)');
ylabel('Im(\Delta t^i_j)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ROOTED SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the complex solutions of the rooted dilation equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(18)
hold on
scatter(real(Rooted_Crosspath_Dilation_Eq_Complex), ...
    imag(Rooted_Crosspath_Dilation_Eq_Complex), 10, 'filled');
scatter(real(Rooted_Crosspath_Dilation_Eq_Complex_Inv), ...
    imag(Rooted_Crosspath_Dilation_Eq_Complex_Inv), 10, 'filled');
title('Rooted dilation equation, complex solutions');
xlabel('Real');
ylabel('Imaginary');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the real only solutions of the rooted dilation equation  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(19)
hold on
scatter(1:size(Rooted_Crosspath_Dilation_Eq_Real, 1),...
    Rooted_Crosspath_Dilation_Eq_Real, 5, 'filled');
scatter(1:size(Rooted_Crosspath_Dilation_Eq_Real_Inv, 1),...
    Rooted_Crosspath_Dilation_Eq_Real_Inv, 5, 'filled');
title('Rooted dilation equation, real only solutions');
xlabel('Pairwise comparison');
ylabel('Real');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Rates of time elapse %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(20)
hold on
scatter(1:size(Rooted_Crosspath_Dilation_Eq_Real, 1),...
    Rooted_Crosspath_Dilation_Eq_Real, 5, 'filled');
title('Rates of time elapse');
xlabel('Pairwise comparison');
ylabel('Relative rate of time elapse');