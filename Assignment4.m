clc;
clear;
close all;

% Toroidal JKR Model: Load vs. Approach for different Rt values

% Parameters
mu_t = 5; % Fixed toroidal Tabor parameter
Rt_values = [100, 50, 20, 10]; % Major radius values for the torus
b_t = linspace(-8, 8, 2000); % Non-dimensional distance bt

% Preallocate matrices for non-dimensional load and approach
w_t_bar = zeros(length(Rt_values), length(b_t));
delta_t_bar = zeros(length(Rt_values), length(b_t));

% Calculate non-dimensional load (w_t_bar) and approach (delta_t_bar) for each Rt
for i = 1:length(Rt_values)
    Rt = Rt_values(i);
    
    % Non-dimensional load (w_t_bar) based on the Toroidal JKR model
    w_t_bar(i, :) = (pi / 4) * (b_t ./ sqrt(mu_t)).^2 - sqrt(2) * sqrt(pi) * (b_t ./ sqrt(mu_t)).^(1/2);
    
    % Non-dimensional approach (delta_t_bar) based on the Toroidal JKR model
    delta_t_bar(i, :) = (1/4) * (b_t ./ sqrt(mu_t)).^2 .* (1 + 2 * log(Rt ./ b_t) + 8 * log(2)) ...
                        - sqrt(2 / pi) * (b_t ./ sqrt(mu_t)).^(1/2) .* (2 * log(Rt ./ b_t) + 8 * log(2));
end

% Plot 1: Load vs. Approach
figure(1);
hold on;
colors = {'r-', 'g--', 'b-', 'm--'};
labels = {
    '$\mu_t = 5, \ \bar{R}_t = 100$', ...
    '$\mu_t = 5, \ \bar{R}_t = 50$', ...
    '$\mu_t = 5, \ \bar{R}_t = 20$', ...
    '$\mu_t = 5, \ \bar{R}_t = 10$'};

% Plot Load vs. Approach for each Rt value
for i = 1:length(Rt_values)
    plot(delta_t_bar(i, :), w_t_bar(i, :), colors{i}, 'LineWidth', 1.5);
end

% Axis labels, title, and legend
legend(labels, 'Interpreter', 'latex', 'Location', 'best');
xlabel('Approach $\frac{\hat{\delta}}{\mu_t} = \frac{\delta}{\mu_t \epsilon}$', 'Interpreter', 'latex');
ylabel('Load $\hat{w} = \frac{w}{2 \pi \Delta \gamma R_t \sqrt{\frac{\epsilon}{\rho}}}$', 'Interpreter', 'latex');
title('Load vs. Approach for Toroidal JKR Model');
grid on;
xlim([-8 2]);
ylim([-2 0.5]);
hold off;

%% Plot 2: Load vs. Contact Half-width
% Preallocation for non-dimensional contact half-width (a_bar)
a_bar = zeros(length(Rt_values), length(w_t_bar));
a_t = b_t ./ sqrt(mu_t);  % Contact half-width

% Plotting Contact Half-width vs. Load (Figure 2)
figure(2);
hold on;

% Plot each Rt value
for i = 1:length(Rt_values)
    plot(w_t_bar(i, :), a_t, colors{i}, 'LineWidth', 1.5);
end

% Axis labels, title, and legend
legend(labels, 'Interpreter', 'latex', 'Location', 'best');
xlabel('Load $\hat{w} = \frac{w}{2 \pi \Delta \gamma R_t \sqrt{\frac{\epsilon}{\rho}}}$', 'Interpreter', 'latex');
ylabel('Contact half-width $\frac{\hat{a}}{\sqrt{\mu_t}} = \frac{a}{\sqrt{\mu_t \epsilon \rho}}$', 'Interpreter', 'latex');
title('Contact Half-width vs. Load for Toroidal JKR Model');
grid on;
xlim([-2 1]);
ylim([0 2.5]);
hold off;

%% Assignment 4.2 Contact splitting and pull off force

% Define constants (example values)
Rt = 10;              % Radius of original torus (e.g., 10 cm)
rho = 1;              % Characteristic length (e.g., 1 cm)
E_star = 1e6;         % Effective Young's modulus (Pa)
delta_gamma = 0.05;   % Surface energy change (J/m^2)

% Function to compute the pull-off force for a toroidal indenter
pull_off_force = @(Rt) 3 * pi * Rt * ((pi * rho * E_star * delta_gamma^2 / 2)^(1/3));

% Logarithmic scale for the number of contacts (from 1 to 1e12)
n_contacts = logspace(0, 12, 100);  % 100 points between 10^0 and 10^12

% Calculate the total pull-off force for each n_contacts
total_pull_off = zeros(size(n_contacts));

for i = 1:length(n_contacts)
    % For each case, split the original contact into 4 smaller contacts
    Rt_new = Rt / 2;  % New radius for each smaller torus (half of original Rt)
    
    % Compute the pull-off force for the smaller toruses
    pull_off_single_contact = pull_off_force(Rt_new);
    
    % Each smaller torus contributes to the total force, scaled by the number of contacts
    total_pull_off(i) = 4 * n_contacts(i) * pull_off_single_contact;  % Multiply by 4 due to the scaling law
end

% Plot the total pull-off force vs number of contacts on a logarithmic scale
figure;

% Plot with a line (without dots) and adjust line style/width
plot(n_contacts, total_pull_off, '-', 'LineWidth', 2, 'Color', [0 0.5 0]);  % Dark green line

% Make the plot a log-log plot
set(gca, 'XScale', 'log', 'YScale', 'log');

% Customize axes and labels
xlabel('Number of contacts (n)');
ylabel('Total pull-off force (N)');
title('Total Adhesion Force vs Number of Contacts for Toroidal JKR Model');

% Customize grid: fewer grid lines and adjust the grid style
grid on;
ax = gca;
ax.GridLineStyle = '--'; % Dashed grid lines
ax.GridAlpha = 0.5;      % Grid transparency
ax.MinorGridLineStyle = ':';  % Minor gridlines are dotted

% Remove the minor gridlines to make the plot cleaner
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

% Adjust plot limits for better clarity (optional)
xlim([1, 1e12]);
ylim([1e-2, max(total_pull_off)]);

% Show the plot
hold off;



