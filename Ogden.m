clear; clc;

% Set the range of lambda values
lambda_min = 0.05;
lambda_max = 2;
[lambda1, lambda2] = meshgrid(linspace(lambda_min, lambda_max, 500));

% Set up a tiled layout with two columns and four rows
figure;
t = tiledlayout(4, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Add row labels for each matrix form
annotation('textbox', [0.03, 0.87, 0.05, 0.05], 'String', 'a)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.03, 0.63, 0.05, 0.05], 'String', 'b)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.03, 0.39, 0.05, 0.05], 'String', 'c)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.03, 0.15, 0.05, 0.05], 'String', 'd)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);

% Define matrix equations based on the images
matrix_definitions = {
    @(lam1, lam2, a1, a2) [
        (lam1^(a1-1) - lam1^(-a1/2-1)), (lam1^(a2-1) - lam1^(-a2/2-1));
        (lam2^(a1-1) - lam2^(-a1/2-1)), (lam2^(a2-1) - lam2^(-a2/2-1))
    ],
    @(lam1, lam2, a1, a2) [
        (lam1^(a1-1) - lam1^(-(a1+1))), (lam1^(a2-1) - lam1^(-a2+1));
        (1 - lam1^(-a1)), (1 - lam1^(-a2));
        (lam2^(a1-1) - lam2^(-(a1+1))), (lam2^(a2-1) - lam2^(-a2+1));
        (1 - lam2^(-a1)), (1 - lam2^(-a2));
    ],
    @(lam1, lam2, a1, a2) [
        (lam1^(a1-1) - lam1^(-2*a1-1)), (lam1^(a2-1) - lam1^(-2*a2-1));
        (lam2^(a1-1) - lam2^(-2*a1-1)), (lam2^(a2-1) - lam2^(-2*a2-1))
    ],
    @(lam1, lam2, a1, a2) [
        (lam1^(a1-1) - lam1^(-a1-1) * lam2^(-a1)), (lam1^(a2-1) - lam1^(-a2-1) * lam2^-(a2));
        (lam2^(a1-1) - lam2^(-a1-1) * lam1^(-a1)), (lam2^(a2-1) - lam2^(-a2-1) * lam1^(-a2))
    ]
};

% Define alpha values
a1_values = [1, 2, 3];
a2_values = [-1, -2, -3];
num_pairs = length(a1_values) * length(a2_values);

% Generate enough colors and line styles for all pairs
contour_colors = lines(num_pairs);
line_styles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-'};

% Set up a tiled layout with two columns and four rows
figure;
t = tiledlayout(4, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Initialize a legend handle array
legend_handles = [];
legend_entries = {};

% Loop through each matrix definition
for matrix_index = 1:length(matrix_definitions)
    matrix_func = matrix_definitions{matrix_index};
    
    pair_index = 1;
    condition_number_sum = zeros(size(lambda1));
    determinant_sum = zeros(size(lambda1));

    % Loop over each alpha pair to compute and plot
    for i = 1:length(a1_values)
        a1 = a1_values(i);
        for j = 1:length(a2_values)
            a2 = a2_values(j);

            % Calculate condition number and determinant for each point in the grid
            condition_number = zeros(size(lambda1));
            determinant = zeros(size(lambda1));
            for m = 1:size(lambda1, 1)
                for n = 1:size(lambda1, 2)
                    lam1 = lambda1(m, n);
                    lam2 = lambda2(m, n);
                    
                    % Define the matrix
                    matrix = matrix_func(lam1, lam2, a1, a2);
                    
                    % Calculate condition number and determinant
                    condition_number(m, n) = sqrt(cond(transpose(matrix) * matrix));
                    determinant(m, n) = sqrt(abs(det(transpose(matrix) * matrix)));
                end
            end

            % Store legend entry
            if matrix_index == 1
                legend_entries{pair_index} = sprintf('\\alpha_1 = %d, \\alpha_2 = %d', a1, a2);
            end
            
            % Accumulate condition numbers and determinants
            condition_number_sum = condition_number_sum + condition_number;
            determinant_sum = determinant_sum + determinant;
            
            % Get current color and style
            current_color = contour_colors(pair_index, :);
            current_style = line_styles{pair_index};
            
            % Plot condition number heatmap and contour
            nexttile(2 * matrix_index - 1);
            hold on;
            avg_condition_number = condition_number_sum / num_pairs;
            h1 = imagesc(lambda1(1, :), lambda2(:, 1), avg_condition_number, 'AlphaData', 0.1);
            set(gca, 'YDir', 'normal');
            colormap(gca, winter);
            clim([1, 20]); % Set color limit for condition number heatmap
            [C1, h_contour1] = contour(lambda1, lambda2, condition_number, [20, 20], 'LineColor', current_color, ...
                                       'LineStyle', current_style, 'LineWidth', 1.5);
            
            % Collect handles and entries for the legend
            if matrix_index == 1
                legend_handles = [legend_handles; h_contour1];
            end
            axis equal; % Make subplot square
            colorbar('Location', 'eastoutside'); % Add colorbar to the side
            box on;
            
            % Plot determinant heatmap and contour
            nexttile(2 * matrix_index);
            hold on;
            avg_determinant = determinant_sum / num_pairs;
            h2 = imagesc(lambda1(1, :), lambda2(:, 1), avg_determinant, 'AlphaData', 0.1);
            set(gca, 'YDir', 'normal');
            colormap(gca, winter);
            clim([0, 1]); % Set color limit for determinant heatmap
            [C2, h_contour2] = contour(lambda1, lambda2, determinant, [1, 1], 'LineColor', current_color, ...
                                       'LineStyle', current_style, 'LineWidth', 1.5);
            axis equal; % Make subplot square
            colorbar('Location', 'eastoutside'); % Add colorbar to the side
            box on;

            pair_index = pair_index + 1;
        end
    end
end

% Add a legend to the first condition number subplot
% nexttile(1);
% legend(legend_handles, legend_entries, 'Location', 'northeast', 'FontSize', 10);
% title('Condition Number Contours Legend');
