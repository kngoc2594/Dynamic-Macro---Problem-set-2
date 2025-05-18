%% File Info.
%{
    simulate.m
    -------
    This code solves the life cycle model with a dynamic budget constraint for two models:
        1. Baseline model: Age-income only
        2. Extended model: Age-icome-education
%}

%% Simulate class.
classdef simulate
    methods(Static)

        %% Baseline Model: Age-Income Only
        function sim = lc(par, sol)
            % Set the random seed for reproducibility
            if isfield(par, 'seed')
                rng(par.seed);
            else
                rng(42);  % Default seed if not specified
                fprintf('Warning: "seed" not specified. Using default value: 42\n');
            end

            %% Parameters and containers for baseline model
            T = par.T;
            NN = par.NN;
            tr = par.tr;
            r = par.r;
            kappa = par.kappa;
            Gt = par.Gt;
            rho = par.rho;
            sigma_eps = par.sigma_eps;
            agrid = par.agrid;

            % Containers for simulation results
            ysim = nan(T, NN);
            asim = nan(T+1, NN);
            csim = nan(T, NN);
            usim = nan(T, NN);

            % Initial income with baseline age effect
            eps0 = randn(1, NN) * sigma_eps;
            ysim(1, :) = Gt(1) .* exp(eps0);
            asim(1, :) = 0;

            % Consumption and saving in the first period
            apol = sol.a;
            cpol = sol.c;
            for i = 1:NN
                age = 1;
                csim(1, i) = cpol(1, age);
                asim(2, i) = apol(1, age);
            end
            usim(1, :) = model.utility(csim(1, :), par);

            % Simulation over time for baseline model
            for t = 2:T
                for i = 1:NN
                    age = t - 1;
                    if age < tr
                        eps_t = randn * sigma_eps;
                        ysim(t, i) = Gt(t) * exp(rho * log(ysim(t-1, i) / Gt(t-1)) + eps_t);
                    else
                        % After retirement: fixed pension based on last working income
                        ysim(t, i) = kappa * ysim(tr-1, i);
                    end

                    % Asset and consumption calculation
                    at = asim(t, i);
                    [~, at_idx] = min(abs(agrid - at));
                    csim(t, i) = cpol(at_idx, t);
                    asim(t+1, i) = apol(at_idx, t);
                    usim(t, i) = model.utility(csim(t, i), par);
                end
            end

            % Final condition: assets at the end of life
            asim(T+1, :) = 0;

            % Output structure for baseline model
            sim = struct();
            sim.ysim = ysim;
            sim.asim = asim(1:T, :);
            sim.csim = csim;
            sim.usim = usim;

            fprintf('\nBaseline Model Simulation Completed.\n');
        end
%% Education-Enhanced Model: Age-Education-Income
function sim = lc_edu(par, sol)
    % Set the random seed for reproducibility
    if isfield(par, 'seed')
        rng(par.seed);
    else
        rng(42);
        fprintf('Warning: "seed" not specified. Using default value: 42\n');
    end

    %% Parameters and containers for education-enhanced model
    T = par.T;
    NN = par.NN;
    tr = par.tr;
    r = par.r;
    kappa = par.kappa;
    Gt = par.Gt;
    rho = par.rho;
    sigma_eps = par.sigma_eps;
    alpha = par.alpha;

    % Ensure the asset grid is generated
    if ~isfield(par, 'agrid')
        par = model.gen_grids(par);
    end
    agrid = par.agrid;

    % Validate Gt length (61 for ages 0 to 60)
    if length(Gt) ~= 61
        error('Error: The Gt vector must have 61 elements (ages 0 to 60).');
    end

    % Check if the education distribution is valid, if not, set default
    if ~isfield(par, 'Edist') || length(par.Edist) ~= 13
        fprintf('Warning: "Edist" not specified or incorrect length. Using uniform distribution.\n');
        par.Edist = ones(1, 13) / 13;
    end

    % Assign education levels and calculate psi
    E = randsample(0:12, NN, true, par.Edist);
    psi = 1 + (E / 12) * alpha;

    % Containers for simulation results
    ysim = nan(T, NN);
    asim = nan(T+1, NN);
    csim = nan(T, NN);
    usim = nan(T, NN);

    % Initial income with education effect
    eps0 = randn(1, NN) * sigma_eps;
    ysim(1, :) = psi .* Gt(1) .* exp(eps0);
    asim(1, :) = 0;

    % Consumption and saving in the first period
    apol = sol.a;
    cpol = sol.c;
    for i = 1:NN
        age = 1;
        csim(1, i) = cpol(1, age);
        asim(2, i) = apol(1, age);
    end
    usim(1, :) = model.utility(csim(1, :), par);

    % Simulation over time for education-enhanced model
    for t = 2:T
        for i = 1:NN
            age = t - 1;
            if age < tr
                eps_t = randn * sigma_eps;
                psi_i = 1 + (E(i) / 12) * alpha;
                % Ensure Gt(t) does not exceed length
                if t > length(Gt)
                    yt = Gt(end);  % Use the last available value if exceeded
                else
                    yt = Gt(t);  % Use the interpolated value
                end
                ysim(t, i) = psi_i * yt * exp(rho * log(ysim(t-1, i) / (psi_i * Gt(t-1))) + eps_t);
            else
                ysim(t, i) = kappa * ysim(tr-1, i);
            end

            % Asset and consumption calculation
            at = asim(t, i);
            [~, at_idx] = min(abs(agrid - at));
            csim(t, i) = cpol(at_idx, t);
            asim(t+1, i) = apol(at_idx, t);
            usim(t, i) = model.utility(csim(t, i), par);
        end
    end

    % Final condition: assets at the end of life
    asim(T+1, :) = 0;

    % Output structure for education-enhanced model
    sim = struct();
    sim.ysim = ysim;
    sim.asim = asim(1:T, :);
    sim.csim = csim;
    sim.usim = usim;
    sim.E = E;
    sim.psi = psi;

    fprintf('\nEducation-Enhanced Model Simulation Completed.\n');
end
end
end
