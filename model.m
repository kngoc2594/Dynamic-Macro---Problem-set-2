%% File Info.
%{
    model.m
    -------
    This code sets up the life-cycle model with a dynamic budget constraint
    driven by Gt and a new education factor Edist
%}

%% Model class.
classdef model
    methods(Static)

        function par = setup()
            par = struct();

            %% Preferences
            par.T = 61;              % Life span: age 0 to T-1 = 60
            par.tr = 41;             % Retirement begins at age 41
            par.beta = 0.94;         % Discount factor (given)
            par.sigma = 2;         % CRRA coefficient (gamma, given)

            %% Prices and Income
            par.r = 0.03;            % Interest rate
            par.kappa = 0.6;         % Pension replacement rate

            %% Income Process Parameters
            par.rho = 0.9;           % Persistence of income shocks
            par.sigma_eps = 0.2;     % Std dev of shock ε_t

            %% Placeholder Gt 
            par.Gt = ones(par.T, 1); % From the empirical Gt data

            %% Simulation Setup
            par.seed = 2025;         % For reproducibility
            par.TT = par.T;          % Simulate full life span
            par.NN = 3000;           % Number of households 
            %% Education Parameters (NEW)
            par.alpha = 0.5;        % Education return scaling
            par.Edist = [];         % Will be loaded from data
        end

        function par = gen_grids(par)
            par.alen = 300;
            par.amax = 30.0;
            par.amin = 0.0;
            par.agrid = linspace(par.amin, par.amax, par.alen)';
        end

        function u = utility(c, par)
            if par.sigma == 1
                u = log(c);
            else
                u = (c.^(1 - par.sigma)) ./ (1 - par.sigma);
            end
        end

    end
end

