%% File Info.
%{
    prepare_income.m
    -------
    This code import the actual data for Income, Age, Gender, Education,
    Wealth, and Consumption and perform data filtering and computation for
    baseline and extended model.
%}
%%
function [Gt_bin, Gt, Gt_matrix, Edist, psi] = prepare_income(filename, alpha)
    % Set default value for alpha 
    if nargin < 2
        alpha = 0.5;  % Default scaling factor
        fprintf('Alpha not provided, using default value: %.2f\n', alpha);
    end
    % Import data
    opts = detectImportOptions(filename);
    data = readtable(filename, opts);
    data.Properties.VariableNames = {'Income', 'Age', 'Gender', 'Education', 'Wealth', 'Consumption'};

    % Fill missing income with 0 (unemployed)
    data.Income(isnan(data.Income)) = 0;

    % Filter only male individuals 
    data = data(data.Gender == 1, :);

    % Create 5-year age bins
    data.AgeBin = floor(data.Age / 5) * 5;

    % Compute log(income); NaN if income = 0
    data.LogIncome = NaN(height(data), 1);
    idx = data.Income > 0;
    data.LogIncome(idx) = log(data.Income(idx));

    % Compute education distribution
    Edist = histcounts(data.Education, -0.5:1:12.5);
    if sum(Edist) == 0
        error('No valid education data found.');
    end
    Edist = Edist / sum(Edist);  % Normalize to get a probability distribution
    
    % Print education distribution for verification
    fprintf('Education Distribution:\n');
    for e = 0:12
        fprintf('Education Level %d: %.2f%%\n', e, Edist(e+1) * 100);
    end

    % Calculate the education multiplier (ψ)
    education_levels = 0:12;
    psi = 1 + (education_levels / 12) * alpha;
    fprintf('\nEducation Multiplier (ψ) for α = %.2f:\n', alpha);
    for e = 0:12
        fprintf('Education Level %d: ψ = %.3f\n', e, psi(e+1));
    end

    % Compute Gt_matrix: age × education specific income
    age_bins = unique(data.AgeBin);
    if length(age_bins) < 2
        error('Insufficient age bins for interpolation. Please check the data.');
    end
    Gt_matrix = zeros(length(age_bins), length(education_levels));

    for a = 1:length(age_bins)
        for e = 0:12
            % Extract the data for the current age and education level
            idx = (data.AgeBin == age_bins(a)) & (data.Education == e);
            logs = data.LogIncome(idx);
            valid_logs = logs(~isnan(logs));  % Remove NaN values
            if ~isempty(valid_logs)
                Gt_matrix(a, e+1) = exp(mean(valid_logs));
            else
                Gt_matrix(a, e+1) = NaN;  % Assign NaN if no valid income data
            end
        end
    end

    % Replace NaNs in Gt_matrix with the overall mean income
    Gt_matrix(isnan(Gt_matrix)) = exp(mean(data.LogIncome(data.LogIncome > -inf)));

    % Compute the baseline Gt (age-specific) as the average across all education levels
    Gt = mean(Gt_matrix, 2, 'omitnan');

    % Create the output table for baseline Gt
    Gt_bin = table(age_bins, Gt, 'VariableNames', {'AgeBin', 'Gt'});

    % Display Gt matrix for verification
    fprintf('\nAge and Education Specific Income (Gt_matrix):\n');
    disp(Gt_matrix);

    fprintf('\nBaseline Age-Specific Income (Gt):\n');
    disp(Gt_bin);
end
