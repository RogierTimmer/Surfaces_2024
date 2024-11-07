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

% Maximum pressure at the center
% p_max = 3 * P / (2 * pi * a^2);  % Hertzian max pressure

% Loop over each ring to compute pressure
for i = 1:n
    if r(i) <= a
        % Hertzian pressure distribution for r <= a
        %p(i) = p_max * (1 - (r(i)^2 / a^2));
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
%analytical deformation outside of contact area according to eq 2.8
b = a - r;
w = zeros(1, n);
for i = 1:n
    if b(i) == 0
        w(i) = ((2*(a^2) - r(i)^2) / (2 * R));
    elseif r(i) > a
        w(i) = (2 * a / (pi * R)) * sqrt(a^2 - b(i)^2) * asin(a / r(i)) * (-1 / (pi * R)) * ((r(i)^2 - b(i)^2) * asin(sqrt(a^2 - b(i)^2) / sqrt(r(i)^2 - b(i)^2)) - (a^2 - b(i)^2));
    else 
        w(i) = 0;
    end
end 
%w = (2 * a / (pi * R)) * sqrt(a^2 - b^2) * asin(a / r) * (-1 / (pi * R)) * ((r^2 - b^2) * asin(sqrt(a^2 - b^2) / sqrt(r^2 - b^2)) - (a^2 - b^2));


%% Pressure distribution for numerical solution
p4 = zeros(1, n);  % Initialize pressure vector
fun = @(x, b, r) (2 * x^2 - b^2) / (sqrt(x^2 - b^2) * sqrt(x^2 - r^2));

% Loop over each ring to compute pressure
for i = 1:n
    if r(i) <= b(i)
        % Pressure distribution for r <= b
        p4(i) = -El / (pi * R) * integral(@(x) fun(x, b(i), r(i)), b(i), a);
        
    elseif b(i) < r(i) && r(i) <= a
        % Pressure distribution for b < r <= a
        p4(i) = -El / (pi * R) * integral(@(x) fun(x, b(i), r(i)), r(i), a);
          
    else
        % Pressure outside contact is 0
        p4(i) = 0;
    end
end
