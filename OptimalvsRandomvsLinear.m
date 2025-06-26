clear; clc;

% Define parameters for the Ogden model
N = 2; % Order of the Ogden model
C = [620.5*2, -689.4*2]; % Coefficients C_p for the Ogden model
alpha = [2, -2]; % Exponents alpha_p for the Ogden model

% Define stretch ratio range
lambda_vals = linspace(1, 2, 64); % Stretch ratios

% Noise levels to iterate over
noise_stds = [0.5, 1, 2, 5]; % Various noise standard deviations

% Number of iterations (random seeds)
num_iterations = 10000;

% Define number of points for random and linear methods
num_points_random = 3;
num_points_linear = 3;

% Initialize storage for results
all_C1_optimal = zeros(num_iterations, length(noise_stds));
all_C2_optimal = zeros(num_iterations, length(noise_stds));
all_C1_random = zeros(num_iterations, length(noise_stds));
all_C2_random = zeros(num_iterations, length(noise_stds));
all_C1_linear = zeros(num_iterations, length(noise_stds));
all_C2_linear = zeros(num_iterations, length(noise_stds));

for n = 1:length(noise_stds)
    noise_std = noise_stds(n); % Current noise level

    for iter = 1:num_iterations
        % Set random seed
        rng(iter);

        % Initialize stress arrays
        uniaxial_stress = zeros(size(lambda_vals));
        uniaxial_stress_gt = zeros(size(lambda_vals));

        % Generate noisy synthetic data for Ogden model
        for i = 1:length(lambda_vals)
            lambda = lambda_vals(i);
            uniaxial_sum = 0;
            for p = 1:N
                uniaxial_sum = uniaxial_sum + C(p) * (lambda^(alpha(p) - 1) - lambda^(-alpha(p)/2 - 1));
            end
            % Add Gaussian noise
            uniaxial_stress(i) = uniaxial_sum + noise_std/100 * uniaxial_sum * randn();
            uniaxial_stress_gt(i) = uniaxial_sum;
        end

        %% Fit Mooney-Rivlin model to UT
        mooney_rivlin_ut = @(lambda, C1, C2) ...
            2 * C1 .* (lambda - 1 ./ lambda.^2) - 2 * C2 .* (1 ./ lambda.^3 - 1);

        % Optimal points for UT
        [~, idx_near_2_ut] = min(abs(lambda_vals - 2));  
        [~, idx_near_1_29_ut] = mink(abs(lambda_vals - 1.29), 2); 
        optimal_indices_ut = [idx_near_2_ut, idx_near_1_29_ut];
        lambda_optimal_ut = lambda_vals(optimal_indices_ut);
        stress_optimal_ut = uniaxial_stress(optimal_indices_ut);

        A_optimal_ut = [2 * (lambda_optimal_ut - 1 ./ lambda_optimal_ut.^2)', ...
                        -2 * (1 ./ lambda_optimal_ut.^3 - 1)'];
        P_optimal_ut = stress_optimal_ut';
        C_optimal_ut = A_optimal_ut \ P_optimal_ut;
        all_C1_optimal(iter, n) = C_optimal_ut(1);
        all_C2_optimal(iter, n) = C_optimal_ut(2);

        % Fit with random points for UT
        random_indices_ut = randperm(length(lambda_vals), num_points_random); 
        lambda_random_ut = lambda_vals(random_indices_ut);
        stress_random_ut = uniaxial_stress(random_indices_ut);

        A_random_ut = [2 * (lambda_random_ut - 1 ./ lambda_random_ut.^2)', ...
                       -2 * (1 ./ lambda_random_ut.^3 - 1)'];
        P_random_ut = stress_random_ut';
        C_random_ut = A_random_ut \ P_random_ut;
        all_C1_random(iter, n) = C_random_ut(1);
        all_C2_random(iter, n) = C_random_ut(2);

        % Fit with linearly spaced points for UT
        linear_indices = round(linspace(1, length(lambda_vals), num_points_linear));
        lambda_linear = lambda_vals(linear_indices);
        stress_linear = uniaxial_stress(linear_indices);

        A_linear = [2 * (lambda_linear - 1 ./ lambda_linear.^2)', ...
                    -2 * (1 ./ lambda_linear.^3 - 1)'];
        P_linear = stress_linear';
        C_linear = A_linear \ P_linear;
        all_C1_linear(iter, n) = C_linear(1);
        all_C2_linear(iter, n) = C_linear(2);
    end
end

% Create one large figure with a tiled layout
fig = figure('Name', 'Distributions of C1 and C2 for Different Noise Levels', 'NumberTitle', 'off');
tiledlayout(length(noise_stds), 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Update legend labels to include number of points
fitting_methods = {['Optimal (3 pts)'], ...
                   ['Random (' num2str(num_points_random) ' pts)'], ...
                   ['Linear (' num2str(num_points_linear) ' pts)']};
line_styles = {'-','--','-.'};
colors = lines(3); % Use MATLAB's built-in color order

for n = 1:length(noise_stds)
    % Extract data for this noise level
    data_C1 = {all_C1_optimal(:, n), all_C1_random(:, n), all_C1_linear(:, n)};
    data_C2 = {all_C2_optimal(:, n), all_C2_random(:, n), all_C2_linear(:, n)};

    %% C1 subplot
    nexttile; 
    hold on;
    h = gobjects(numel(data_C1),1); % store plot handles for legend
    for i = 1:numel(data_C1)
        [f, xi] = ksdensity(data_C1{i});
        h(i) = plot(xi, f, 'Color', colors(i,:), 'LineStyle', line_styles{i}, 'LineWidth', 1);
    end
    title(['C1 Dist. (Noise STD = ', num2str(noise_stds(n)), '% )']);
    xlabel('C1');
    ylabel('Density');
    if n == length(noise_stds)
        % Only show legend on the first subplot (top-left)
        legend(h, fitting_methods, 'Location', 'best');
    end
    hold off;

    % Limit x-axis for C1 subplot
    all_data_C1 = vertcat(data_C1{:});
    mean_C1 = mean(all_data_C1);
    std_C1 = std(all_data_C1);
    xlim([mean_C1 - 3*std_C1, mean_C1 + 3*std_C1]);
    ylim([0, 30*10^-3]); % Limit y-axis

    %% C2 subplot
    nexttile; 
    hold on;
    for i = 1:numel(data_C2)
        [f, xi] = ksdensity(data_C2{i});
        plot(xi, f, 'Color', colors(i,:), 'LineStyle', line_styles{i}, 'LineWidth', 1);
    end
    title(['C2 Dist. (Noise STD = ', num2str(noise_stds(n)), '% )']);
    xlabel('C2');
    ylabel('Density');
    hold off;

    % Limit x-axis for C2 subplot
    all_data_C2 = vertcat(data_C2{:});
    mean_C2 = mean(all_data_C2);
    std_C2 = std(all_data_C2);
    xlim([mean_C2 - 3*std_C2, mean_C2 + 3*std_C2]);
    ylim([0, 30*10^-3]); % Limit y-axis
end
