% Set the range of lambda values
lambda_min = 0.05;
lambda_max = 2;

% Define a grid for lambda1 and lambda2 for the first three pairs
[lambda1, lambda2] = meshgrid(linspace(lambda_min, lambda_max, 200));

% Set up a tiled layout with two columns and four rows, with custom padding
figure;
t = tiledlayout(4, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Add row labels with more spacing from the plots
annotation('textbox', [0.03, 0.87, 0.05, 0.05], 'String', 'a)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.03, 0.63, 0.05, 0.05], 'String', 'b)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.03, 0.39, 0.05, 0.05], 'String', 'c)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.03, 0.15, 0.05, 0.05], 'String', 'd)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);

% First Pair
compute_cond1 = @(lam1, lam2) sqrt(cond(transpose([2 * (lam1 - 1 / lam1^2);
                                                   2 * (lam2 - 1 / lam2^2)])*...
                                                  [2 * (lam1 - 1 / lam1^2);
                                                   2 * (lam2 - 1 / lam2^2)]));
det1 = @(lam1, lam2) sqrt(abs(det(transpose([2 * (lam1 - 1 / lam1^2);
                                                   2 * (lam2 - 1 / lam2^2)])*...
                                                  [2 * (lam1 - 1 / lam1^2);
                                                   2 * (lam2 - 1 / lam2^2)])));
condition_number1 = arrayfun(compute_cond1, lambda1, lambda2);
determinant1 = arrayfun(det1, lambda1, lambda2);

% Condition Number Plot (First Pair)
nexttile;
surf(lambda1, lambda2, condition_number1, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda^1'); ylabel('\lambda^2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([1, 50]);
title('SQRT of Condition Number');

% Determinant Plot (First Pair)
nexttile;
surf(lambda1, lambda2, determinant1, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda^1'); ylabel('\lambda^2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([0, 1]);
title('SQRT of Determinant');

% Second Pair
compute_cond2 = @(lam1, lam2) sqrt(cond(transpose([2 * (lam1 - 1 / lam1^3);
                                                   2 * (1 - 1 / lam1^2);
                                                   2 * (lam2 - 1 / lam2^3); ...
                                                   2 * (1 - 1 / lam2^2)]) * ...
                                                  [2 * (lam1 - 1 / lam1^3);
                                                   2 * (1 - 1 / lam1^2);
                                                   2 * (lam2 - 1 / lam2^3); ...
                                                   2 * (1 - 1 / lam2^2)]));
det2 = @(lam1, lam2) sqrt(abs(det(transpose([2 * (lam1 - 1 / lam1^3);
                                             2 * (1 - 1 / lam1^2);
                                             2 * (lam2 - 1 / lam2^3); ...
                                             2 * (1 - 1 / lam2^2)]) * ...
                                            [2 * (lam1 - 1 / lam1^3);
                                             2 * (1 - 1 / lam1^2);
                                             2 * (lam2 - 1 / lam2^3); ...
                                             2 * (1 - 1 / lam2^2)])));
condition_number2 = arrayfun(compute_cond2, lambda1, lambda2);
determinant2 = arrayfun(det2, lambda1, lambda2);

% Condition Number Plot (Second Pair)
nexttile;
surf(lambda1, lambda2, condition_number2, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda^1'); ylabel('\lambda^2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([1, 50]);

% Determinant Plot (Second Pair)
nexttile;
surf(lambda1, lambda2, determinant2, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda^1'); ylabel('\lambda^2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([0, 1]);

% Third Pair
compute_cond3 = @(lam1, lam2) sqrt(cond(transpose([2 * (lam1 - 1 / lam1^5);
                                                   2 * (lam2 - 1 / lam2^5)]) * ...
                                                  [2 * (lam1 - 1 / lam1^5);
                                                   2 * (lam2 - 1 / lam2^5)]));
det3 = @(lam1, lam2) sqrt(abs(det(transpose([2 * (lam1 - 1 / lam1^5);
                                             2 * (lam2 - 1 / lam2^5)]) * ...
                                            [2 * (lam1 - 1 / lam1^5);
                                             2 * (lam2 - 1 / lam2^5)])));
condition_number3 = arrayfun(compute_cond3, lambda1, lambda2);
determinant3 = arrayfun(det3, lambda1, lambda2);

% Condition Number Plot (Third Pair)
nexttile;
surf(lambda1, lambda2, condition_number3, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda^1'); ylabel('\lambda^2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([1, 50]);

% Determinant Plot (Third Pair)
nexttile;
surf(lambda1, lambda2, determinant3, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda^1'); ylabel('\lambda^2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([0, 1]);

% Fourth Pair
compute_cond4 = @(lam1, lam2) sqrt(cond(transpose([2 * (lam1 - 1 / (lam1^3 * lam2));
                                                   2 * (lam2 - 1 / (lam1 * lam2^3))]) * ...
                                                  [2 * (lam1 - 1 / (lam1^3 * lam2));
                                                   2 * (lam2 - 1 / (lam1 * lam2^3))]));
det4 = @(lam1, lam2) sqrt(abs(det(transpose([2 * (lam1 - 1 / (lam1^3 * lam2));
                                             2 * (lam2 - 1 / (lam1 * lam2^3))]) * ...
                                            [2 * (lam1 - 1 / (lam1^3 * lam2));
                                             2 * (lam2 - 1 / (lam1 * lam2^3))])));
condition_number4 = arrayfun(compute_cond4, lambda1, lambda2);
determinant4 = arrayfun(det4, lambda1, lambda2);

% Condition Number Plot (Third Pair)
nexttile;
surf(lambda1, lambda2, condition_number4, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda_1'); ylabel('\lambda_2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([1, 50]);

% Determinant Plot (Third Pair)
nexttile;
surf(lambda1, lambda2, determinant4, 'EdgeColor', 'none');
view(2); colormap jet; colorbar; axis square;
xlabel('\lambda_1'); ylabel('\lambda_2'); xlim([lambda_min, lambda_max]); ylim([lambda_min, lambda_max]); clim([0, 1]);