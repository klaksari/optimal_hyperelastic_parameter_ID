clear; clc;

% Define parameters
num_points = 64; % Number of points for lambda in each dimension
lambda_values = linspace(0.75, 2, num_points); % Generate lambda values from 0.4 to 1.0

% Maximum number of tests for each loading mode
max_n_uniaxial = 4;       % Uniaxial tests
max_n_equibiaxial = 2;    % Equibiaxial tests
max_n_biaxial = 1;        % Biaxial tests
max_n_pure_shear = 2;     % Pure shear tests

%% Initialize data structures for each loading mode
uniaxial_results = {};
equibiaxial_results = {};
biaxial_results = {};
pure_shear_results = {};

%% Function to pad arrays
pad = @(x, n) [x; NaN(n - numel(x), 1)];

%% Uniaxial Test Analysis
for n = 1:max_n_uniaxial
    % Generate all combinations of lambda values for n tests
    lambda_grid = cell(1, n);
    [lambda_grid{:}] = ndgrid(lambda_values);
    lambda_combinations = cell2mat(cellfun(@(x) x(:), lambda_grid, 'UniformOutput', false));

    % Initialize counters and arrays
    num_combinations = size(lambda_combinations, 1);
    count_cond_gt_20 = 0;
    count_det_lt_1 = 0;
    count_both_conditions = 0;
    det_values = zeros(num_combinations, 1);
    cond_values = zeros(num_combinations, 1);

    % Loop through each combination of lambda values
    for i = 1:num_combinations
        lambdas = lambda_combinations(i, :);

        % Mooney-Rivlin Uniaxial Test Analysis
        A_uniaxial_mr = zeros(n, 2);
        for j = 1:n
            lam1 = lambdas(j);
            A_uniaxial_mr(j, 1) = 2 * (lam1 - 1 / lam1^2);
            A_uniaxial_mr(j, 2) = -2 * (1 / lam1^3 - 1);
        end
        ATA_uniaxial_mr = A_uniaxial_mr' * A_uniaxial_mr;
        det_uniaxial_mr = sqrt(abs(det(ATA_uniaxial_mr)));
        cond_uniaxial_mr = sqrt(cond(ATA_uniaxial_mr));

        % Store det and cond values
        det_values(i) = det_uniaxial_mr;
        cond_values(i) = cond_uniaxial_mr;

        % Count configurations where cond > 20
        if cond_uniaxial_mr > 20
            count_cond_gt_20 = count_cond_gt_20 + 1;
        end

        % Count configurations where det < 1
        if det_uniaxial_mr < 1
            count_det_lt_1 = count_det_lt_1 + 1;
        end

        % Count configurations where both conditions are true
        if cond_uniaxial_mr > 20 && det_uniaxial_mr < 1
            count_both_conditions = count_both_conditions + 1;
        end
    end

    % Calculate the percentages
    percentage_cond_gt_20 = (count_cond_gt_20 / num_combinations) * 100;
    percentage_det_lt_1 = (count_det_lt_1 / num_combinations) * 100;
    percentage_both_conditions = (count_both_conditions / num_combinations) * 100;

    % Calculate statistical values
    max_det = max(det_values);
    min_det = min(det_values);
    avg_det = mean(det_values);
    std_det = std(det_values);
    median_det = median(det_values);

    max_cond = max(cond_values);
    min_cond = min(cond_values);
    avg_cond = mean(cond_values);
    std_cond = std(cond_values);
    median_cond = median(cond_values);

    % Find lambda combinations
    min_cond_indices = find(cond_values == min_cond);
    lambda_min_cond = lambda_combinations(min_cond_indices(1), :);

    max_det_indices = find(det_values == max_det);
    lambda_max_det = lambda_combinations(max_det_indices(1), :);

    % Calculate the expression for each combination
    expr_values = (cond_values / median_cond).^2 + (median_det ./ det_values).^2;
    expr_values(det_values == 0) = inf; % Avoid division by zero

    % Find the lambda combination minimizing the expression
    min_expr = min(expr_values);
    min_expr_indices = find(expr_values == min_expr);
    lambda_min_expr = lambda_combinations(min_expr_indices(1), :);

    % Retrieve det and cond at lambda_min_expr
    det_min_expr = det_values(min_expr_indices(1));
    cond_min_expr = cond_values(min_expr_indices(1));

    % Convert lambda combinations to strings
    lambda_min_cond_str = mat2str(lambda_min_cond);
    lambda_max_det_str = mat2str(lambda_max_det);
    lambda_min_expr_str = mat2str(lambda_min_expr);

    % Store results
    uniaxial_results = [uniaxial_results; {n, percentage_cond_gt_20, percentage_det_lt_1, percentage_both_conditions, ...
                        max_det, min_det, avg_det, std_det, median_det, lambda_max_det_str, ...
                        max_cond, min_cond, avg_cond, std_cond, median_cond, lambda_min_cond_str, ...
                        lambda_min_expr_str, det_min_expr, cond_min_expr}];
end

%% Equibiaxial Test Analysis
for n = 1:max_n_equibiaxial
    % Generate all combinations of lambda values for n tests
    lambda_grid = cell(1, n);
    [lambda_grid{:}] = ndgrid(lambda_values);
    lambda_combinations = cell2mat(cellfun(@(x) x(:), lambda_grid, 'UniformOutput', false));

    % Initialize counters and arrays
    num_combinations = size(lambda_combinations, 1);
    count_cond_gt_50 = 0;
    count_det_lt_1 = 0;
    count_both_conditions = 0;
    det_values = zeros(num_combinations, 1);
    cond_values = zeros(num_combinations, 1);

    % Loop through each combination of lambda values
    for i = 1:num_combinations
        lambdas = lambda_combinations(i, :);

        % Mooney-Rivlin Equibiaxial Test Analysis
        A_equibiaxial_mr = zeros(n, 2);
        for j = 1:n
            lam1 = lambdas(j);
            A_equibiaxial_mr(j, 1) = 2 * (lam1 - 1 / lam1^5);
            A_equibiaxial_mr(j, 2) = -2 * (1 / lam1^3 - lam1^3);
        end
        ATA_equibiaxial_mr = A_equibiaxial_mr' * A_equibiaxial_mr;
        det_equibiaxial_mr = sqrt(abs(det(ATA_equibiaxial_mr)));
        cond_equibiaxial_mr = sqrt(cond(ATA_equibiaxial_mr));

        % Store det and cond values
        det_values(i) = det_equibiaxial_mr;
        cond_values(i) = cond_equibiaxial_mr;

        % Count configurations where cond > 20
        if cond_equibiaxial_mr > 20
            count_cond_gt_20 = count_cond_gt_20 + 1;
        end

        % Count configurations where det < 1
        if det_equibiaxial_mr < 1
            count_det_lt_1 = count_det_lt_1 + 1;
        end

        % Count configurations where both conditions are true
        if cond_equibiaxial_mr > 20 && det_equibiaxial_mr < 1
            count_both_conditions = count_both_conditions + 1;
        end
    end

    % Calculate statistical values
    percentage_cond_gt_20 = (count_cond_gt_20 / num_combinations) * 100;
    percentage_det_lt_1 = (count_det_lt_1 / num_combinations) * 100;
    percentage_both_conditions = (count_both_conditions / num_combinations) * 100;

    max_det = max(det_values);
    min_det = min(det_values);
    avg_det = mean(det_values);
    std_det = std(det_values);
    median_det = median(det_values);

    max_cond = max(cond_values);
    min_cond = min(cond_values);
    avg_cond = mean(cond_values);
    std_cond = std(cond_values);
    median_cond = median(cond_values);

    % Find lambda combinations
    min_cond_indices = find(cond_values == min_cond);
    lambda_min_cond = lambda_combinations(min_cond_indices(1), :);

    max_det_indices = find(det_values == max_det);
    lambda_max_det = lambda_combinations(max_det_indices(1), :);

    % Calculate the expression
    expr_values = (cond_values / median_cond) + (median_det ./ det_values);
    expr_values(det_values == 0) = inf;

    % Find lambda minimizing the expression
    min_expr = min(expr_values);
    min_expr_indices = find(expr_values == min_expr);
    lambda_min_expr = lambda_combinations(min_expr_indices(1), :);

    % Retrieve det and cond at lambda_min_expr
    det_min_expr = det_values(min_expr_indices(1));
    cond_min_expr = cond_values(min_expr_indices(1));

    % Convert to strings
    lambda_min_cond_str = mat2str(lambda_min_cond);
    lambda_max_det_str = mat2str(lambda_max_det);
    lambda_min_expr_str = mat2str(lambda_min_expr);

    % Store results
    equibiaxial_results = [equibiaxial_results; {n, percentage_cond_gt_20, percentage_det_lt_1, percentage_both_conditions, ...
                             max_det, min_det, avg_det, std_det, median_det, lambda_max_det_str, ...
                             max_cond, min_cond, avg_cond, std_cond, median_cond, lambda_min_cond_str, ...
                             lambda_min_expr_str, det_min_expr, cond_min_expr}];
end

%% Biaxial Test Analysis
for n = 1:max_n_biaxial
    % Generate all combinations of lambda_1 and lambda_2 values for n tests
    lambda_grid = cell(1, 2 * n);
    [lambda_grid{:}] = ndgrid(lambda_values);
    lambda_combinations = cell2mat(cellfun(@(x) x(:), lambda_grid, 'UniformOutput', false));

    % Initialize counters and arrays
    num_combinations = size(lambda_combinations, 1);
    count_cond_gt_50 = 0;
    count_det_lt_1 = 0;
    count_both_conditions = 0;
    det_values = zeros(num_combinations, 1);
    cond_values = zeros(num_combinations, 1);

    % Loop through each combination of lambda values
    for i = 1:num_combinations
        lambdas = reshape(lambda_combinations(i, :), [2, n]);

        % Mooney-Rivlin Biaxial Test Analysis
        A_biaxial_mr = zeros(2 * n, 2);
        for j = 1:n
            lam1 = lambdas(1, j);
            lam2 = lambdas(2, j);
            A_biaxial_mr(2*j-1, 1) = 2 * (lam1 - 1 / (lam1^3 * lam2^2));
            A_biaxial_mr(2*j-1, 2) = -2 * (1 / lam1^3 - lam1 * lam2^2);
            A_biaxial_mr(2*j,   1) = 2 * (lam2 - 1 / (lam1^2 * lam2^3));
            A_biaxial_mr(2*j,   2) = -2 * (1 / lam2^3 - lam1^2 * lam2);
        end
        ATA_biaxial_mr = A_biaxial_mr' * A_biaxial_mr;
        det_biaxial_mr = sqrt(abs(det(ATA_biaxial_mr)));
        cond_biaxial_mr = sqrt(cond(ATA_biaxial_mr));

        % Store det and cond values
        det_values(i) = det_biaxial_mr;
        cond_values(i) = cond_biaxial_mr;

        % Count configurations
        if cond_biaxial_mr > 20
            count_cond_gt_20 = count_cond_gt_20 + 1;
        end
        if det_biaxial_mr < 1
            count_det_lt_1 = count_det_lt_1 + 1;
        end
        if cond_biaxial_mr > 20 && det_biaxial_mr < 1
            count_both_conditions = count_both_conditions + 1;
        end
    end

    % Calculate statistical values
    percentage_cond_gt_20 = (count_cond_gt_20 / num_combinations) * 100;
    percentage_det_lt_1 = (count_det_lt_1 / num_combinations) * 100;
    percentage_both_conditions = (count_both_conditions / num_combinations) * 100;

    max_det = max(det_values);
    min_det = min(det_values);
    avg_det = mean(det_values);
    std_det = std(det_values);
    median_det = median(det_values);

    max_cond = max(cond_values);
    min_cond = min(cond_values);
    avg_cond = mean(cond_values);
    std_cond = std(cond_values);
    median_cond = median(cond_values);

    % Find lambda combinations
    min_cond_indices = find(cond_values == min_cond);
    lambda_min_cond = lambda_combinations(min_cond_indices(1), :);

    max_det_indices = find(det_values == max_det);
    lambda_max_det = lambda_combinations(max_det_indices(1), :);

    % Calculate the expression
    expr_values = (cond_values / median_cond) + (median_det ./ det_values);
    expr_values(det_values == 0) = inf;

    % Find lambda minimizing the expression
    min_expr = min(expr_values);
    min_expr_indices = find(expr_values == min_expr);
    lambda_min_expr = lambda_combinations(min_expr_indices(1), :);

    % Retrieve det and cond at lambda_min_expr
    det_min_expr = det_values(min_expr_indices(1));
    cond_min_expr = cond_values(min_expr_indices(1));

    % Convert to strings
    lambda_min_cond_str = mat2str(lambda_min_cond);
    lambda_max_det_str = mat2str(lambda_max_det);
    lambda_min_expr_str = mat2str(lambda_min_expr);

    % Store results
    biaxial_results = [biaxial_results; {n, percentage_cond_gt_20, percentage_det_lt_1, percentage_both_conditions, ...
                          max_det, min_det, avg_det, std_det, median_det, lambda_max_det_str, ...
                          max_cond, min_cond, avg_cond, std_cond, median_cond, lambda_min_cond_str, ...
                          lambda_min_expr_str, det_min_expr, cond_min_expr}];
end

%% Pure Shear Test Analysis
for n = 1:max_n_pure_shear
    % Generate all combinations of lambda values for n tests
    lambda_grid = cell(1, n);
    [lambda_grid{:}] = ndgrid(lambda_values);
    lambda_combinations = cell2mat(cellfun(@(x) x(:), lambda_grid, 'UniformOutput', false));

    % Initialize counters and arrays
    num_combinations = size(lambda_combinations, 1);
    count_cond_gt_50 = 0;
    count_det_lt_1 = 0;
    count_both_conditions = 0;
    det_values = zeros(num_combinations, 1);
    cond_values = zeros(num_combinations, 1);

    % Loop through all combinations of lambda values
    for i = 1:num_combinations
        lambdas = lambda_combinations(i, :);

        % Initialize the matrix A for this combination
        A_pure_shear_mr = zeros(2 * n, 2);
        for j = 1:n
            lam = lambdas(j);

            % Construct the rows of A based on the equations
            A_pure_shear_mr(2 * j - 1, 1) = 2 * (lam - 1 / lam^3);             
            A_pure_shear_mr(2 * j - 1, 2) = -2 * (1 / lam^3 - lam);           
            A_pure_shear_mr(2 * j, 1) = 2 * (1 - 1 / lam^2);                  
            A_pure_shear_mr(2 * j, 2) = -2 * (1 - lam^2);                  
        end

        % Calculate ATA, determinant, and condition number
        ATA_pure_shear_mr = A_pure_shear_mr' * A_pure_shear_mr;
        det_pure_shear_mr = sqrt(abs(det(ATA_pure_shear_mr)));
        cond_pure_shear_mr = sqrt(cond(ATA_pure_shear_mr));

        % Store determinant and condition number values
        det_values(i) = det_pure_shear_mr;
        cond_values(i) = cond_pure_shear_mr;

        % Count configurations
        if cond_pure_shear_mr > 20
            count_cond_gt_20 = count_cond_gt_20 + 1;
        end
        if det_pure_shear_mr < 1
            count_det_lt_1 = count_det_lt_1 + 1;
        end
        if cond_pure_shear_mr > 20 && det_pure_shear_mr < 1
            count_both_conditions = count_both_conditions + 1;
        end
    end

    % Calculate statistical values
    percentage_cond_gt_20 = (count_cond_gt_20 / num_combinations) * 100;
    percentage_det_lt_1 = (count_det_lt_1 / num_combinations) * 100;
    percentage_both_conditions = (count_both_conditions / num_combinations) * 100;

    max_det = max(det_values);
    min_det = min(det_values);
    avg_det = mean(det_values);
    std_det = std(det_values);
    median_det = median(det_values);

    max_cond = max(cond_values);
    min_cond = min(cond_values);
    avg_cond = mean(cond_values);
    std_cond = std(cond_values);
    median_cond = median(cond_values);

    % Find lambda combinations for specific conditions
    min_cond_indices = find(cond_values == min_cond);
    lambda_min_cond = lambda_combinations(min_cond_indices(1), :);

    max_det_indices = find(det_values == max_det);
    lambda_max_det = lambda_combinations(max_det_indices(1), :);

    % Calculate the expression
    expr_values = (cond_values / median_cond) + (median_det ./ det_values);
    expr_values(det_values == 0) = inf;

    % Find lambda minimizing the expression
    min_expr = min(expr_values);
    min_expr_indices = find(expr_values == min_expr);
    lambda_min_expr = lambda_combinations(min_expr_indices(1), :);

    % Retrieve det and cond at lambda_min_expr
    det_min_expr = det_values(min_expr_indices(1));
    cond_min_expr = cond_values(min_expr_indices(1));

    % Convert to strings for storage
    lambda_min_cond_str = mat2str(lambda_min_cond);
    lambda_max_det_str = mat2str(lambda_max_det);
    lambda_min_expr_str = mat2str(lambda_min_expr);

    % Store results
    pure_shear_results = [pure_shear_results; {n, percentage_cond_gt_20, percentage_det_lt_1, percentage_both_conditions, ...
                             max_det, min_det, avg_det, std_det, median_det, lambda_max_det_str, ...
                             max_cond, min_cond, avg_cond, std_cond, median_cond, lambda_min_cond_str, ...
                             lambda_min_expr_str, det_min_expr, cond_min_expr}];
end


%% Create tables for each loading mode
var_names = {'Number_of_Tests_n', 'Cond_GT_20_%', 'Det_LT_1_%', 'Both_Conditions_%', ...
             'Max_Det', 'Min_Det', 'Avg_Det', 'Std_Det', 'Median_Det', 'Lambda_Max_Det', ...
             'Max_Cond', 'Min_Cond', 'Avg_Cond', 'Std_Cond', 'Median_Cond', 'Lambda_Min_Cond', ...
             'Lambda_Min_Expr', 'Det_at_Lambda_Min_Expr', 'Cond_at_Lambda_Min_Expr'};

uniaxial_table = cell2table(uniaxial_results, 'VariableNames', var_names);
equibiaxial_table = cell2table(equibiaxial_results, 'VariableNames', var_names);
biaxial_table = cell2table(biaxial_results, 'VariableNames', var_names);
pure_shear_table = cell2table(pure_shear_results, 'VariableNames', var_names);

% Display the tables
disp('Uniaxial Test Results (Mooney-Rivlin Model)');
disp(uniaxial_table);

disp('Equibiaxial Test Results (Mooney-Rivlin Model)');
disp(equibiaxial_table);

disp('Biaxial Test Results (Mooney-Rivlin Model)');
disp(biaxial_table);

disp('Pure Shear Test Results (Mooney-Rivlin Model)');
disp(pure_shear_table);

%% Save the tables to CSV files
min_lambda = min(lambda_values);
max_lambda = max(lambda_values);

output_folder = sprintf('Mooney_Rivlin_tables_lambda_%.2f_to_%.2f', min_lambda, max_lambda);
output_folder = strrep(output_folder, '.', '_');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save results
writetable(uniaxial_table, fullfile(output_folder, 'uniaxial_results.csv'));
writetable(equibiaxial_table, fullfile(output_folder, 'equibiaxial_results.csv'));
writetable(biaxial_table, fullfile(output_folder, 'biaxial_results.csv'));
writetable(pure_shear_table, fullfile(output_folder, 'pure_shear_results.csv'));

%% Plotting Results with Subplots in One Figure

% Extract data for plotting from the tables
% For each loading mode, extract n, max_det, median_det, min_cond, median_cond

% Uniaxial Test Data
n_uniaxial = [uniaxial_table.Number_of_Tests_n];
max_det_uniaxial = [uniaxial_table.Max_Det];
median_det_uniaxial = [uniaxial_table.Median_Det];
min_cond_uniaxial = [uniaxial_table.Min_Cond];
median_cond_uniaxial = [uniaxial_table.Median_Cond];

% Equibiaxial Test Data
n_equibiaxial = [equibiaxial_table.Number_of_Tests_n];
max_det_equibiaxial = [equibiaxial_table.Max_Det];
median_det_equibiaxial = [equibiaxial_table.Median_Det];
min_cond_equibiaxial = [equibiaxial_table.Min_Cond];
median_cond_equibiaxial = [equibiaxial_table.Median_Cond];

% Biaxial Test Data
n_biaxial = [biaxial_table.Number_of_Tests_n];
max_det_biaxial = [biaxial_table.Max_Det];
median_det_biaxial = [biaxial_table.Median_Det];
min_cond_biaxial = [biaxial_table.Min_Cond];
median_cond_biaxial = [biaxial_table.Median_Cond];

% Pure Shear Test Data
n_pure_shear = [pure_shear_table.Number_of_Tests_n];
max_det_pure_shear = [pure_shear_table.Max_Det];
median_det_pure_shear = [pure_shear_table.Median_Det];
min_cond_pure_shear = [pure_shear_table.Min_Cond];
median_cond_pure_shear = [pure_shear_table.Median_Cond];

%% Combined Figure: Determinants and Condition Numbers
figure('Name', 'Determinants and Condition Numbers');
max_n_overall = max([max_n_uniaxial, max_n_equibiaxial, ...
                     max_n_biaxial,  max_n_pure_shear]);
% Adjust the figure size (optional)
set(gcf, 'Position', [100, 100, 1200, 600]); % [left, bottom, width, height]

% Define color order for consistency
colors = lines(4); % 4 different colors for the 4 tests

%% Subplot 1: Determinants (Median & Maximum)
subplot(1, 2, 1); % 1 row, 2 columns, first plot
xticks(1:max_n_overall);
hold on;

% Plotting for Uniaxial Test
plot(n_uniaxial, median_det_uniaxial, 'o--', 'Color', colors(1, :), 'DisplayName', 'Uniaxial (Median)', 'LineWidth', 1.5);
plot(n_uniaxial, max_det_uniaxial, 'o-', 'Color', colors(1, :), 'DisplayName', 'Uniaxial (Max)', 'LineWidth', 1.5);

% Plotting for Equibiaxial Test
plot(n_equibiaxial, median_det_equibiaxial, 's--', 'Color', colors(2, :), 'DisplayName', 'Equibiaxial (Median)', 'LineWidth', 1.5);
plot(n_equibiaxial, max_det_equibiaxial, 's-', 'Color', colors(2, :), 'DisplayName', 'Equibiaxial (Max)', 'LineWidth', 1.5);

% Plotting for Biaxial Test
plot(n_biaxial, median_det_biaxial, 'd--', 'Color', colors(3, :), 'DisplayName', 'Biaxial (Median)', 'LineWidth', 1.5);
plot(n_biaxial, max_det_biaxial, 'd-', 'Color', colors(3, :), 'DisplayName', 'Biaxial (Max)', 'LineWidth', 1.5);

% Plotting for Pure Shear Test
plot(n_pure_shear, median_det_pure_shear, '^--', 'Color', colors(4, :), 'DisplayName', 'Pure Shear (Median)', 'LineWidth', 1.5);
plot(n_pure_shear, max_det_pure_shear, '^-', 'Color', colors(4, :), 'DisplayName', 'Pure Shear (Max)', 'LineWidth', 1.5);

xlabel('# of Measurements (n)');
ylabel('SQRT Determinant Value');
title('Determinants: Median (Dashed) and Max (Solid)');
%legend('Location', 'best');
grid on;
ylim([0 50]);
hold off;

%% Subplot 2: Condition Numbers (Median & Minimum)
subplot(1, 2, 2); % 1 row, 2 columns, second plot
xticks(1:max_n_overall);
hold on;

% Plotting for Uniaxial Test
plot(n_uniaxial, median_cond_uniaxial, 'o--', 'Color', colors(1, :), 'DisplayName', 'Uniaxial (Median)', 'LineWidth', 1.5);
plot(n_uniaxial, min_cond_uniaxial, 'o-', 'Color', colors(1, :), 'DisplayName', 'Uniaxial (Min)', 'LineWidth', 1.5);

% Plotting for Equibiaxial Test
plot(n_equibiaxial, median_cond_equibiaxial, 's--', 'Color', colors(2, :), 'DisplayName', 'Equibiaxial (Median)', 'LineWidth', 1.5);
plot(n_equibiaxial, min_cond_equibiaxial, 's-', 'Color', colors(2, :), 'DisplayName', 'Equibiaxial (Min)', 'LineWidth', 1.5);

% Plotting for Biaxial Test
plot(n_biaxial, median_cond_biaxial, 'd--', 'Color', colors(3, :), 'DisplayName', 'Biaxial (Median)', 'LineWidth', 1.5);
plot(n_biaxial, min_cond_biaxial, 'd-', 'Color', colors(3, :), 'DisplayName', 'Biaxial (Min)', 'LineWidth', 1.5);

% Plotting for Pure Shear Test
plot(n_pure_shear, median_cond_pure_shear, '^--', 'Color', colors(4, :), 'DisplayName', 'Pure Shear (Median)', 'LineWidth', 1.5);
plot(n_pure_shear, min_cond_pure_shear, '^-', 'Color', colors(4, :), 'DisplayName', 'Pure Shear (Min)', 'LineWidth', 1.5);

xlabel('# of Measurements (n)');
ylabel('SQRT Condition Number');
title('Condition Numbers: Median (Dashed) and Min (Solid)');
%legend('Location', 'best');
grid on;
ylim([1 50]);
hold off;
