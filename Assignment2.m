clc
clear
close all

% Parameters for the problem
E1 = 200e9;
v1 = 0.3;
E2 = 3e9;
v2 = 0.33;
P = 10;  % Load (N)
El = (E1/(1-v1^2))+(E2/(1-v2^2));  % Reduced modulus of elasticity (Pa)

% Variables for numerical solution
R = 0.05;  % Reduced radius (m)
a = ((3*P*R)/(4*El))^(1/3);  % Contact radius (m)

n = 500;  % Number of rings
r = linspace(0, 2*a, n);  % Radial positions

%% Pressure distribution for numerical solution
p = zeros(1, n);  % Initialize pressure vector

% Loop over each ring to compute pressure
for i = 1:n
    if r(i) <= a
        % Hertzian pressure distribution for r <= a
        p(i)= (2*El)/(pi*R)*sqrt(a^2 - r(i)^2);
    else
        % Pressure outside contact area is zero
        p(i) = 0;
    end
end

%% Assignment 2.2 Compute numerical displacement using the function
K = makeK_constantpressuresphericalring(r, El);
uz_numerical = K * p';

% Analytical Hertzian solution for uz(r)
uz_analytical = zeros(1, n);
r_analytical = linspace(0, 2*a, n);  % Radial positions
for i = 1:n
    if r_analytical(i) <= a
        uz_analytical(i) = ((2*(a^2) - r_analytical(i)^2) / (2 * R));
    else
        uz_analytical(i) = 0;
    end
end 

% Plot numerical vs analytical results
figure;
plot(r, uz_numerical, 'r', 'LineWidth', 2);
hold on
plot(r, uz_analytical, 'b', 'LineWidth', 2);
legend('Numerical Solution',"Analytical Solution");
xlabel('Radius r (m)');
ylabel('Displacement uz (m)');
title('Hertzian Displacement');
grid on;
%% Assignment 2.3, finding the pressure vector from a known displacement and K-matrix
p3 = inv(K) * uz_analytical';

%plotting the difference between the pressure in 2.2  and 2.3
figure;
plot(r, p, 'r', 'LineWidth', 2);
hold on
plot(r,p3,"b","LineWidth",2)
legend('P-situation 2',"P-situation 3");
xlabel('Radius r (m)');
ylabel('Pressure (Pa)');
title('Pressure distribution');
grid on;

%% Assignment 2.4
% MATLAB code for comparing analytical and numerical deformation due to JKR pressure profile

% Parameters
b = 5e-5;          % Radius for singularity region (m)
num_points = 250;  % Number of points in the radial direction

% Generate radial points inside and outside the contact area
r_inside = linspace(0, a, num_points);      % r <= a
r_outside = linspace(a + 1e-6, 2*a, num_points); % r > a

% Analytical pressure profile inside contact area (Eq. 2.79)
p_analytical = zeros(1, num_points);
for i = 1:num_points
    r = r_inside(i);
    if r <= b
        % Integrate the function using numerical integration
        integrand = @(x) (2 * x.^2 - b^2) ./ (sqrt(x.^2 - b^2) .* sqrt(x.^2 - r^2));
        p_analytical(i) = -El / (pi * R) * integral(integrand, b, a);
    elseif r > b && r <= a
        integrand = @(x) (2 * x.^2 - b^2) ./ (sqrt(x.^2 - b^2) .* sqrt(x.^2 - r^2));
        p_analytical(i) = -El / (pi * R) * integral(integrand, r, a);
    end
end

% Analytical deformation inside and outside contact area
w_inside = zeros(1, num_points);
for i = 1:num_points
    r = r_inside(i);
    w_inside(i) = uz_analytical(i); % analytical deformation calculated in 2.2 for inside deformation
end

w_outside = zeros(1, num_points);
for i = 1:num_points
    r = r_outside(i);
    term1 = (2 * a / (pi * R)) * sqrt(a^2 - b^2) * asin(a / r);
    term2 = (1 / (pi * R)) * ((r^2 - b^2) * asin(sqrt((a^2 - b^2) / (r^2 - b^2))) ...
        - sqrt(a^2 - b^2) * sqrt(r^2 - a^2));
    w_outside(i) = term1 - term2; % analytical deformation calculated from eq. 2.79 and 2.80
end

% Combine radial and deformation arrays for a single plot
r_combined = [r_inside, r_outside];
w_analytical_combined = [w_inside, w_outside];

% Plotting the combined deformation results
figure;
plot(r_combined, w_analytical_combined, 'b-', 'LineWidth', 1.5);
hold on;
plot(r_combined, uz_numerical, 'g--', 'LineWidth', 1.5);
xlabel('Radius (r)');
ylabel('Deformation w(r)');
title('Combined Analytical and Numerical Deformation');
legend('Analytical', 'Numerical');
grid on;

%% Assignment 2.5
% Parameters
a = linspace(0, 100, 5000); % Normalized contact radius values
delta_gamma_values = [4,7,10]; % Example values for Delta gamma
zeta = 3; % Assuming force-controlled experiment
E_star = 2.2315; % [Pa], Effective modulus, bit random number 
theta = 30; %random theta in degrees

% Convert angle to radians
theta_rad = deg2rad(theta);

% Initialize figures
figure;
hold on;
title('Normalized Force (F/Fc) and Indentation Depth (d/dc) vs. Normalized Contact Radius (a/ac)');
xlabel('a/a_c');
ylabel('Normalized Values');
grid on;
xlim([0 3]);
ylim([-1 7]);

figure;
hold on;
title('Normalized Force (F/Fc) vs. Normalized Indentation Depth (d/dc)');
xlabel('d/d_c');
ylabel('F/F_c');
grid on;
xlim([-1 3]);
ylim([-1 0.5]);

% Loop over each delta_gamma to calculate and plot normalized values
for i = 1:length(delta_gamma_values)
    delta_gamma = delta_gamma_values(i);
    
    % Calculate critical contact radius (a_c) for the current delta_gamma
    a_c = (2 * zeta^2 * delta_gamma) / (pi * E_star * (tan(theta_rad))^2);
    
    % Define normalized contact radius a/a_c for current delta_gamma
    a_normalized = a / a_c;
    
    % Calculate d/d_c (normalized indentation depth) using the provided equation
    d_normalized = 3 * a_normalized - 2 * sqrt(a_normalized);
    
    % Calculate F/F_c (normalized force) using the provided equation
    F_normalized = 3 * a_normalized.^2 - 4 * sqrt(a_normalized.^3);
    
    % Plot F/Fc and d/dc vs a/ac for each delta_gamma in Figure 6
    figure(4);
    plot(a_normalized, F_normalized, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('F/Fc, \\Delta\\gamma = %.1f', delta_gamma));
    plot(a_normalized, d_normalized, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('d/dc, \\Delta\\gamma = %.1f', delta_gamma));
    
    % Plot F/Fc vs d/dc for each delta_gamma in Figure 7
    figure(5);
    plot(d_normalized, F_normalized, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('\\Delta\\gamma = %.1f', delta_gamma));
end

% Finalize the legends for the plots
figure(4);
legend('Location', 'northwest');
hold off;

figure(5);
legend('Location', 'southwest');
hold off;

