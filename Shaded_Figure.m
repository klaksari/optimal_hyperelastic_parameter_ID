function plot_shaded_mr_cross_prediction_with_frequency()
% This script:
%   1) Generates synthetic uniaxial data using the Mooney–Rivlin model with
%      ground‐truth parameters (C1 = 620.5, C2 = 689.4) plus noise.
%   2) Obtains Mooney–Rivlin fits (via three strategies: Optimal, Random, Linear)

%   3) Uses these fits to predict the stress in three deformation modes:
%         Uniaxial, Equibiaxial, and Pure Shear.
%   4) Displays two figures:
%         (a) A 3×3 grid (rows: modes, columns: strategies) where each subplot 
%             shows the individual repeated experiments as a scatter plot along 
%             with the ground‐truth response. The y‐axis is cut off at 5000.
%         (b) A 3×3 grid of shaded regions. In each subplot the envelope (i.e. the
%             area between the minimum and maximum predicted stress curves) is
%             plotted, along with the ground‐truth response. The y‐axis is cut
%             off at 5000 to better show differences.
%
% In all plots the ground‐truth response is computed using the Mooney–Rivlin model
% with parameters C1 = 620.5 and C2 = 689.4.

    % -------------
    % USER SETTINGS
    % -------------
    num_iterations     = 20;  % Number of Monte Carlo repeats
    noise_std_percent  = 5;    % 5% noise in uniaxial data
    num_points_random  = 3;    % # points for random strategy
    num_points_linear  = 3;    % # points for linear strategy

    % Ground truth Mooney–Rivlin parameters
    C1_true = 620.5;
    C2_true = 689.4;

    % The stretch ratios for uniaxial data generation & fitting
    lambda_vals_uniax = linspace(0.75, 2, 64);

    % "Optimal" point choices (for uniaxial test data)
    lambda_opt_target1 = 2;
    lambda_opt_target2 = 0.75;

    %===========================
    % 1) Generate & Fit Uniaxial Data (Monte Carlo)
    %===========================
    [C1_opt_array, C2_opt_array, ...
     C1_rand_array, C2_rand_array, ...
     C1_lin_array,  C2_lin_array] = run_uniaxial_fits( ...
        num_iterations, noise_std_percent, ...
        num_points_random, num_points_linear, ...
        lambda_vals_uniax, ...
        lambda_opt_target1, lambda_opt_target2, ...
        C1_true, C2_true);

    %===========================
    % 2) Define the stretch range for cross-prediction
    %===========================
    lambda_range = linspace(0.75, 2, 100);

    %===========================
    % 3) Build the Cross-Prediction Figure (3 rows x 3 columns)
    %    (Rows: Uniaxial, Equibiaxial, Pure Shear)
    %    (Cols: Optimal (green), Random (red), Linear (blue))
    %    In each subplot the individual Monte Carlo predictions are plotted as scatter points
    %    (with low marker opacity) and the ground‐truth curve is overlaid in dashed magenta.
    %    The y-axis is limited to 5000.
    %===========================
    figure('Name','Cross-Validations from UT Fits','NumberTitle','off');
    set(gcf, 'Position',[100,100,1400,900]); % Optional window sizing

    % --- Row 1: Uniaxial (the fitted mode) ---
    subplot(3,3,1); % Optimal
    plot_scatter_mode(lambda_range, 'uniaxial', ...
        C1_opt_array, C2_opt_array, 'g', 'Optimal (Uniaxial)', C1_true, C2_true);

    subplot(3,3,2); % Random
    plot_scatter_mode(lambda_range, 'uniaxial', ...
        C1_rand_array, C2_rand_array, 'r', 'Random (Uniaxial)', C1_true, C2_true);

    subplot(3,3,3); % Linear
    plot_scatter_mode(lambda_range, 'uniaxial', ...
        C1_lin_array, C2_lin_array, 'b', 'Linear (Uniaxial)', C1_true, C2_true);

    % --- Row 2: Equibiaxial Cross-Prediction ---
    subplot(3,3,4); % Optimal
    plot_scatter_mode(lambda_range, 'equibiaxial', ...
        C1_opt_array, C2_opt_array, 'g', 'Optimal (Equibiaxial Cross - Validation)', C1_true, C2_true);

    subplot(3,3,5); % Random
    plot_scatter_mode(lambda_range, 'equibiaxial', ...
        C1_rand_array, C2_rand_array, 'r', 'Random (Equibiaxial Cross - Validation)', C1_true, C2_true);

    subplot(3,3,6); % Linear
    plot_scatter_mode(lambda_range, 'equibiaxial', ...
        C1_lin_array, C2_lin_array, 'b', 'Linear (Equibiaxial Cross - Validation)', C1_true, C2_true);

    % --- Row 3: Pure Shear Cross-Prediction ---
    subplot(3,3,7); % Optimal
    plot_scatter_mode(lambda_range, 'pure_shear', ...
        C1_opt_array, C2_opt_array, 'g', 'Optimal (Pure Shear Cross - Validation)', C1_true, C2_true);

    subplot(3,3,8); % Random
    plot_scatter_mode(lambda_range, 'pure_shear', ...
        C1_rand_array, C2_rand_array, 'r', 'Random (Pure Shear Cross - Validation)', C1_true, C2_true);

    subplot(3,3,9); % Linear
    plot_scatter_mode(lambda_range, 'pure_shear', ...
        C1_lin_array, C2_lin_array, 'b', 'Linear (Pure Shear Cross - Validation)', C1_true, C2_true);

    %===========================
    % 4) Build the Frequency Map Figure (Shaded Region, 3 rows x 3 columns)
    %    Each subplot shows a shaded region that encloses the stress predictions
    %    (i.e. the envelope formed by the minimum and maximum values at each λ)
    %    with the ground‐truth curve overlaid.
    %    The y-axis is limited to 5000.
    %===========================
    figure('Name','Frequency Map with Shaded Region','NumberTitle','off');
    set(gcf, 'Position',[100,100,1400,900]); % Optional window sizing

    modeNames   = {'uniaxial','equibiaxial','pure_shear'};
    modeTitles  = {'Uniaxial','Equibiaxial Cross-Validation','Pure Shear Cross-Validation'};
    stratLabels = {'Optimal','Random','Linear'};

    % Store the fitted parameter arrays for convenience
    all_C1 = {C1_opt_array, C1_rand_array, C1_lin_array};
    all_C2 = {C2_opt_array, C2_rand_array, C2_lin_array};

    for iMode = 1:3
        modeName  = modeNames{iMode};
        modeTitle = modeTitles{iMode};

        for jStrat = 1:3
            c1Array   = all_C1{jStrat};
            c2Array   = all_C2{jStrat};
            stratLabel = stratLabels{jStrat};

            subplotIndex = (iMode-1)*3 + jStrat;
            subplot(3,3,subplotIndex);
            hold on; grid on; box on;

            numIter   = length(c1Array);
            numLambda = length(lambda_range);
            stressMat = zeros(numIter, numLambda);

            % Compute stress predictions for each Monte Carlo iteration
            for iter = 1:numIter
                stressMat(iter,:) = mooney_rivlin_stress_mode(modeName, lambda_range, c1Array(iter), c2Array(iter));
            end

            % Compute the envelope (minimum and maximum stress at each lambda)
            stressMin = min(stressMat, [], 1);
            stressMax = max(stressMat, [], 1);

            % Create the patch (shaded region) data.
            Xpatch = [lambda_range, fliplr(lambda_range)];
            Ypatch = [stressMin, fliplr(stressMax)];

            % Choose a color based on the strategy.
            switch lower(stratLabel)
                case 'optimal'
                    patchColor = [0, 1, 0]; % green
                case 'random'
                    patchColor = [1, 0, 0]; % red
                case 'linear'
                    patchColor = [0, 0, 1]; % blue
                otherwise
                    patchColor = [0.5, 0.5, 0.5];
            end

            % Plot the shaded region.
            fill(Xpatch, Ypatch, patchColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Overlay the ground‐truth curve.
            gt_stress = mooney_rivlin_stress_mode(modeName, lambda_range, C1_true, C2_true);
            plot(lambda_range, gt_stress, '--m', 'LineWidth',2, 'DisplayName','Ground Truth');

            title(sprintf('%s (%s)', modeTitle, stratLabel), 'FontWeight','bold');
            xlabel('\lambda');
            ylabel('Cauchy Stress');
            ylim([0, 7000]);  % Cut off the y-axis at 5000 for clarity.
            legend('Location','best');
            hold off;
        end
    end
end

% ===================================================================
% Subfunctions
% ===================================================================

function [C1_opt_array, C2_opt_array, ...
          C1_rand_array, C2_rand_array, ...
          C1_lin_array,  C2_lin_array] = run_uniaxial_fits( ...
    num_iterations, noise_std_percent, ...
    num_points_random, num_points_linear, ...
    lambda_vals, lambda_opt_target1, lambda_opt_target2, ...
    C1_true, C2_true)
% Generates uniaxial data (with noise) using the Mooney–Rivlin model with
% ground‐truth parameters (C1_true, C2_true) then fits Mooney–Rivlin parameters
% for each iteration using three strategies: Optimal, Random, and Linear.

    C1_opt_array  = zeros(num_iterations,1);
    C2_opt_array  = zeros(num_iterations,1);
    C1_rand_array = zeros(num_iterations,1);
    C2_rand_array = zeros(num_iterations,1);
    C1_lin_array  = zeros(num_iterations,1);
    C2_lin_array  = zeros(num_iterations,1);

    for iter = 1:num_iterations
        rng(iter);  % For reproducibility

        % 1) Generate synthetic uniaxial stress data using the Mooney–Rivlin model
        uniaxial_stress = zeros(size(lambda_vals));
        for i = 1:length(lambda_vals)
            lam = lambda_vals(i);
            sigma_mr = 2 * C1_true * (lam - 1./lam.^2) - 2 * C2_true * (1./lam.^3 - 1);
            noise = (noise_std_percent/100) * sigma_mr * randn();
            uniaxial_stress(i) = sigma_mr + noise;
        end

        % 2) Define the Mooney–Rivlin uniaxial formula
        mooney_uni = @(l, C1, C2) 2 * C1 .* (l - 1./l.^2) - 2 * C2 .* (1./l.^3 - 1);

        % 3) Fit (A) Optimal: Choose 2 specific points near λ = lambda_opt_target1 and λ = lambda_opt_target2
        [~, idx_near_2]   = min(abs(lambda_vals - lambda_opt_target1));
        [~, idx_near_129] = mink(abs(lambda_vals - lambda_opt_target2), 2);
        idx_optimal       = [idx_near_2, idx_near_129];

        lam_opt = lambda_vals(idx_optimal);
        P_opt   = uniaxial_stress(idx_optimal);

        A_opt = [ 2*(lam_opt - 1./lam_opt.^2)', ...
                 -2*(1./lam_opt.^3 - 1)' ];
        C_opt = A_opt \ P_opt';
        C1_opt_array(iter) = C_opt(1);
        C2_opt_array(iter) = C_opt(2);

        % 4) Fit (B) Random: Pick num_points_random random data points
        rand_idx = randperm(length(lambda_vals), num_points_random);
        lam_rand = lambda_vals(rand_idx);
        P_rand   = uniaxial_stress(rand_idx);

        A_rand = [ 2*(lam_rand - 1./lam_rand.^2)', ...
                  -2*(1./lam_rand.^3 - 1)' ];
        C_rand = A_rand \ P_rand';
        C1_rand_array(iter) = C_rand(1);
        C2_rand_array(iter) = C_rand(2);

        % 5) Fit (C) Linear: Pick num_points_linear equally spaced points,
        %    excluding the first point (lambda = 1)
        linear_indices = round(linspace(2, length(lambda_vals), num_points_linear));
        lam_lin = lambda_vals(linear_indices);
        % Print the selected lambdas for linear fit
        fprintf('Iteration %d: Selected lambdas for linear fit: %s\n', iter, mat2str(lam_lin));
        P_lin   = uniaxial_stress(linear_indices);

        A_lin = [ 2*(lam_lin - 1./lam_lin.^2)', ...
                 -2*(1./lam_lin.^3 - 1)' ];
        C_lin = A_lin \ P_lin';
        C1_lin_array(iter) = C_lin(1);
        C2_lin_array(iter) = C_lin(2);
    end
end

% -------------------------------------------------------------------
function plot_scatter_mode(lambda_vals, modeName, C1_array, C2_array, colorChar, plotLabel, C1_true, C2_true)
% Plots a scatter plot of the repeated stress predictions for a given mode.
% Each individual experiment’s prediction is plotted as a scatter point.
% The ground‐truth curve (computed via the Mooney–Rivlin model) is overlaid.
%
% Inputs:
%   modeName   - 'uniaxial', 'equibiaxial', or 'pure_shear'
%   C1_array   - Array of fitted C1 values from Monte Carlo
%   C2_array   - Array of fitted C2 values from Monte Carlo
%   colorChar  - Color for the scatter markers (e.g., 'g', 'r', 'b')
%   plotLabel  - Title for the subplot
%   C1_true    - Ground‐truth C1 value (620.5)
%   C2_true    - Ground‐truth C2 value (689.4)

    hold on; grid on; box on;

    numIter   = length(C1_array);
    numLambda = length(lambda_vals);
    stressMat = zeros(numIter, numLambda);

    % Compute stress for each Monte Carlo iteration using the fitted parameters
    for i = 1:numIter
        stressMat(i,:) = mooney_rivlin_stress_mode(modeName, lambda_vals, ...
            C1_array(i), C2_array(i));
    end

    % Create X coordinates for scatter: replicate lambda values for each experiment.
    X = repmat(lambda_vals, numIter, 1);
    
    % Plot all stress predictions as scatter points.
    scatter(X(:), stressMat(:), 5, colorChar, 'filled', 'MarkerFaceAlpha', 0.3);
    
    % Overlay the ground‐truth curve computed via the Mooney–Rivlin model
    gt_stress = mooney_rivlin_stress_mode(modeName, lambda_vals, C1_true, C2_true);
    plot(lambda_vals, gt_stress, '--k', 'LineWidth',2, 'DisplayName','Ground Truth');
    
    title(plotLabel, 'FontWeight','bold');
    xlabel('\lambda');
    ylabel('Cauchy Stress');
    
    % switch lower(modeName)
    %     case 'uniaxial'
    %         ylim([0,10000]);
    %     case 'equibiaxial'
    %         ylim([-5000,100000]);
    %     case 'pure_shear'
    %         ylim([0,10000]);
    %     otherwise
    %         error('Unknown modeName: %s', modeName);
    % end
    hold off;
end

% -------------------------------------------------------------------
function stress_vals = mooney_rivlin_stress_mode(modeName, lambda_vals, C1, C2)
% Computes the Mooney–Rivlin **Cauchy** stress for the specified deformation mode.
%
% Uniaxial     : σ = 2 C1 (λ − λ⁻²) + 2 C2 (1 − λ⁻³)
% Equibiaxial  : σ = 2 (C1 + C2 λ²) (λ² − λ⁻⁴)
% Pure shear   : σ = 2 (C1 + C2) (λ² − λ⁻²)
%
% Inputs
%   modeName     'uniaxial' | 'equibiaxial' | 'pure_shear'
%   lambda_vals  vector of stretches λ
%   C1, C2       Mooney–Rivlin material parameters
%
% Output
%   stress_vals  vector of Cauchy stresses (same size as lambda_vals)

    switch lower(modeName)
        case 'uniaxial'           % λ1 = λ,  λ2 = λ3 = λ⁻½
            stress_vals = 2*C1.*(lambda_vals       - lambda_vals.^-2) ...
                        + 2*C2.*(1                - lambda_vals.^-3);

        case 'equibiaxial'        % λ1 = λ2 = λ,  λ3 = λ⁻²
            stress_vals = 2*(C1 + C2.*lambda_vals.^2) ...
                        .* (lambda_vals.^2 - lambda_vals.^-4);

        case 'pure_shear'         % λ1 = λ,  λ2 = λ⁻¹,  λ3 = 1
            stress_vals = 2*(C1 + C2) ...
                        .* (lambda_vals.^2 - lambda_vals.^-2);

        otherwise
            error('Unknown modeName: %s', modeName);
    end
end
