clc
clear
close all

% Assignment 3 - Surfaces - Rogier, Niels and Dirk

%% 3.1

Rq = 0.1e-6; % root-mean-square roughness [0.1 micrometer]

H = 0.5; % Hurst exponent
Lx = 1e-5; % Topography length [10 micrometer]
m = 100; % nr of pixels in x direction
n = 100; % nr of pixels in y direction

% in artificial_surf, nanmean(A) had to be adjusted to
% mean(A,"omitmissing") in order to make it run
[z,PixelWidth, PSD] = artificial_surf(Rq,H,Lx,m,n);


% plotting the surface
x = linspace(0,Lx,m);
y = linspace(0,(Lx*n/m),n);
s = surf(x,y,z);
s.EdgeColor = 'interp';

%% 3.2 a)
% APDF

% plotting the amplitude probability density function
height_values = z(:);
num_bins = 250; % Number of bins, adjust based on resolution
[counts, bin_edges] = histcounts(height_values, num_bins, 'Normalization', 'pdf');

% Calculate the bin centers
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% Plot the PDF
figure;
plot(bin_centers, counts, 'LineWidth', 2, 'DisplayName', 'Empirical PDF');
hold on;

% Calculate Gaussian distribution (Normal Distribution)
mu = mean(height_values);    % Mean of the height values
sigma = std(height_values);  % Standard deviation of the height values

% Generate Gaussian distribution based on mean and standard deviation
gaussian_pdf = (1/(sigma * sqrt(2 * pi))) * exp(-(bin_centers - mu).^2 / (2 * sigma^2));

% Plot the Gaussian distribution
plot(bin_centers, gaussian_pdf, 'r--', 'LineWidth', 2, 'DisplayName', 'Gaussian Distribution');

% Add labels and title
xlabel('Amplitude (Height Values)');
ylabel('Probability Density');
title('Amplitude Probability Density Function with Gaussian Distribution');
legend('show');
grid on;
hold off;

%% 3.2 b)
% ACF

pixel_size = Lx*1e6 / m;  % Size of each pixel in micrometers (total length/pixels in x-direction)

% Select a specific row of height values (ACF in the x-direction)
row_idx = 25;  % Selecting the middle row
height_slice = z(row_idx, :);  % Extracting heights along the x-direction at row 25

% Calculate the autocorrelation function (ACF)
[acf, lags] = xcorr(height_slice, 'coeff');  % Normalized autocorrelation
positive_lags = lags(lags >= 0);             % Select only positive lags
acf_positive = acf(lags >= 0);               % Corresponding ACF values for positive lags

% Convert lags from pixels to micrometers
positive_lags_micrometers = positive_lags * pixel_size;

% Fit an exponential decay to the ACF
% Model: ACF(lag) = exp(-lag / L), where L is the autocorrelation length
fit_func = @(b, x) exp(-x / b);               % Exponential decay function
initial_guess = 5;                            % Initial guess for the autocorrelation length
fit_params = nlinfit(positive_lags_micrometers, acf_positive, fit_func, initial_guess);  % Fit the model

% Calculate the autocorrelation length from the fit
autocorr_length = fit_params;

% Generate the exponential fit line for plotting
exp_fit_curve = fit_func(autocorr_length, positive_lags_micrometers);

% Generate x-axis values for the height slice in micrometers
x_values_micrometers = (0:m-1) * pixel_size;

% Plot the height values of the slice and the ACF with exponential fit side by side
figure;

% Plot 1: Height values of the slice
subplot(1, 2, 1);  % Create the first subplot in a 1x2 grid
plot(x_values_micrometers, height_slice, 'LineWidth', 2);
xlabel('x (μm)');
ylabel('Height');
title('Height Values of the Slice');
grid on;

% Plot 2: ACF with the exponential fit
subplot(1, 2, 2);  % Create the second subplot in a 1x2 grid
plot(positive_lags_micrometers, acf_positive, 'LineWidth', 2, 'DisplayName', 'ACF');
hold on;
plot(positive_lags_micrometers, exp_fit_curve, 'r--', 'LineWidth', 2, ...
    'DisplayName', ['Exponential Fit (L = ' num2str(autocorr_length, '%.2f') ' \mum)']);
xlabel('Lag (μm)');
ylabel('Autocorrelation');
title('Autocorrelation Function (ACF) with Exponential Fit');
legend('show');
grid on;
hold off;

% Set the figure title
sgtitle('Height Values and Autocorrelation Function (ACF) for a slice in x-direction');

%% 3.2 c)
% SF

% Calculate the variance A(0)
A0 = mean(height_slice.^2);

% Initialize the structure function array
S_tau = zeros(1, m);

% Calculate the structure function S(tau) using the autocovariance approach
for tau = 0:m-1
    % Calculate the autocovariance A(tau)
    if tau == 0
        A_tau = A0;
    else
        % Shifted product mean for A(tau)
        A_tau = mean(height_slice(1:end-tau) .* height_slice(1+tau:end));
    end
    
    % Calculate the structure function using the formula
    S_tau(tau + 1) = 2 * (A0 - A_tau);
end

% Define the x-axis (displacement values)
delta_x = (0:m-1) * PixelWidth;  % PixelWidth is the spacing in the x-direction

% Plot the structure function
figure;
plot(delta_x, S_tau, 'LineWidth', 2);
xlabel('\tau [m]')
ylabel('S(\tau)')
title('Structure Function in the x-direction (Autocovariance Method)')
grid on;

%% 3.3 1)
% show fractal behacior

% Define the x-axis (displacement values) and exclude zero to avoid log(0)
delta_x = (0:m-1) * PixelWidth;
nonzero_idx = delta_x > 0;  % Ignore zero lag

% Log-log plot of S(tau) vs. delta_x (to confirm power-law behavior)
log_delta_x = log(delta_x(nonzero_idx));
log_S_tau = log(S_tau(nonzero_idx));

% Plot structure function on a log-log scale
figure;
loglog(delta_x, S_tau, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Structure Function S(\tau)');
hold on;

% Fit a linear model on the log-log data to find the slope
p = polyfit(log_delta_x, log_S_tau, 1);

% Plot the fitted line to visualize the scaling
loglog(delta_x(nonzero_idx), exp(polyval(p, log_delta_x)), '--', 'DisplayName', sprintf('Fit: Slope = %.2f', p(1)));

% Add labels and title
xlabel('\tau [m]');
ylabel('S(\tau)');
title('Structure Function with Log-Log Plot and Power Law Fit');
legend;
grid on;

% Display the slope of the fit, which should be close to 1 if H = 0.5
disp(['Estimated slope (2H) = ', num2str(p(1))]);

%% 3.3 2)
% show that the surface is isotropic

% Structure function comparison

% Select a specific row (x-direction) and column (y-direction) slice
row_idx = round(m/2);  % Selecting the middle row for the x-direction
col_idx = round(n/2);  % Selecting the middle column for the y-direction

% Extract height values along the x-direction and y-direction
height_slice_x = z(row_idx, :);  % Along the x-direction at row 25
height_slice_y = z(:, col_idx);  % Along the y-direction at column 25

% Calculate the variance A(0) for both directions
A0_x = mean(height_slice_x.^2);
A0_y = mean(height_slice_y.^2);

% Initialize the structure function arrays
S_tau_x = zeros(1, m);  % Structure function in x-direction
S_tau_y = zeros(1, n);  % Structure function in y-direction

% Calculate the structure function S(tau) for the x-direction
for tau = 0:m-1
    if tau == 0
        A_tau_x = A0_x;
    else
        A_tau_x = mean(height_slice_x(1:end-tau) .* height_slice_x(1+tau:end));
    end
    S_tau_x(tau + 1) = 2 * (A0_x - A_tau_x);
end

% Calculate the structure function S(tau) for the y-direction
for tau = 0:n-1
    if tau == 0
        A_tau_y = A0_y;
    else
        A_tau_y = mean(height_slice_y(1:end-tau) .* height_slice_y(1+tau:end));
    end
    S_tau_y(tau + 1) = 2 * (A0_y - A_tau_y);
end

% Define the displacement values in both directions
delta_x = (0:m-1) * PixelWidth;
delta_y = (0:n-1) * PixelWidth;

% Plot the structure functions in both directions
figure;
plot(delta_x, S_tau_x, 'r-', 'LineWidth', 2, 'DisplayName', 'S_x(\tau)');
hold on;
plot(delta_y, S_tau_y, 'b-', 'LineWidth', 2, 'DisplayName', 'S_y(\tau)');
xlabel('\tau [m]')
ylabel('S(\tau)')
title('Structure Function in x- and y-directions')
legend;
grid on;

%% 3.2 2)
% Different approach
% Radially averaged PSD

% Compute the 2D Fourier Transform of the surface
Z_fft = fft2(z);
PSD_2D = abs(Z_fft).^2;

% Shift the zero-frequency component to the center
PSD_2D_shifted = fftshift(PSD_2D);

% Define the grid for the radial averaging
[rows, cols] = size(z);
center_row = floor(rows / 2) + 1;
center_col = floor(cols / 2) + 1;
[X, Y] = meshgrid(1:cols, 1:rows);
R = sqrt((X - center_col).^2 + (Y - center_row).^2);

% Radially average the PSD
max_radius = floor(min(rows, cols) / 2);
radial_psd = zeros(1, max_radius);
for r = 1:max_radius
    radial_psd(r) = mean(PSD_2D_shifted(R >= r-0.5 & R < r+0.5), 'all');
end

% Plot the radially averaged PSD
figure;
plot((1:max_radius) * (1 / Lx), radial_psd, 'LineWidth', 1.5);
xlabel('Spatial Frequency [1/m]');
ylabel('Radially Averaged PSD');
title('Radially Averaged Power Spectral Density');
grid on;