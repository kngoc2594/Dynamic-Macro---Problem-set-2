%% File Info.
%{
    solve.m
    -------
    This code solves the life cycle model with a dynamic budget constraint,
    using backward induction and age-specific stochastic income (Gt).
%}

%% Solve class.
classdef solve
    methods(Static)

        function sol = lc(par)
            %% Structure array for model solution
            sol = struct();

            %% Model parameters, grids, and functions
            T = par.T;
            tr = par.tr;
            beta = par.beta;
            kappa = par.kappa;
            r = par.r;
            Gt = par.Gt;

            % Check for 'agrid' in the structure, set default if missing
            if isfield(par, 'agrid')
                agrid = par.agrid;
            else
                % Set a default asset grid if not provided
                agrid = linspace(0, 10, 100);  
                fprintf('Warning: "agrid" not specified. Using default grid of length: %d\n', length(agrid));
            end

            % Set the length of the asset grid
            alen = length(agrid);

            %% Backward induction
            v1 = nan(alen, T);  % Value function
            a1 = nan(alen, T);  % Optimal asset choice
            c1 = nan(alen, T);  % Optimal consumption

            fprintf('------------Solving from the Last Period of Life.------------\n\n')

            % Backward induction: Solve from the last period to the first
            for age = T-1:-1:0
                t = age + 1; 

            if t == T
                % Last period: consume all available resources
                if tr > length(Gt)
                    yt = kappa * Gt(end);  
                    fprintf('Warning: tr (%d) exceeds Gt length (%d). Using Gt(end) for pension calculation.\n', tr, length(Gt));
                else
                    yt = kappa * Gt(tr);   
                end
                c1(:, t) = agrid + yt;
                a1(:, t) = 0.0;
                v1(:, t) = model.utility(c1(:, t), par);
            else
            % Income calculation based on retirement status
            if t >= tr
            % Check for tr does not exceed the length of 'Gt'
                if tr > length(Gt)
                    yt = kappa * Gt(end);  % Use the last available Gt value for pension
                    fprintf('Warning: Retirement age threshold (tr) exceeds Gt length. Using Gt(end) for pension calculation.\n');
                else
                    yt = kappa * Gt(tr);  % Pension income after retirement
                end
            else
                % Ensure 't' does not exceed the length of 'Gt'
                if t > length(Gt)
                    yt = Gt(end);  % Use the last available Gt value for income
                    fprintf('Warning: t (%d) exceeds Gt length (%d). Using Gt(end) for income calculation.\n', t, length(Gt));
                else
                    yt = Gt(t);  % Working income
                end
            end
            % Loop over each asset level in the grid
                for i = 1:alen
                at = agrid(i);  % Current asset level
                % Compute possible consumption given assets and income
                ct = at + yt - (agrid ./ (1 + r));
                ct(ct <= 0) = NaN;  % Eliminate non-positive consumption

                % Calculate utility for feasible consumption
                util = model.utility(ct, par);
                next_v = v1(:, t+1);  % Future value function

                % Bellman equation: current utility + discounted future value
                val = util + beta * next_v;
                val(isnan(ct)) = -inf;  % Exclude infeasible values

                % Ensure 'val' is a column vector
                if isvector(val)
                    [vmax, ind] = max(val);  % Find maximum value and index
                else
                    [vmax, ind] = max(val(:));  % Flatten to ensure a scalar result
                end

                % Assign the scalar maximum value
                v1(i, t) = vmax;
                c1(i, t) = ct(ind);
                a1(i, t) = agrid(ind);  % Optimal asset choice
                    end
                end

                % Progress update every 5 periods
                if mod(t, 5) == 0
                    fprintf('Solved Age: %d\n', age);
                end
            end

            %% Output the solution
            sol.c = c1;
            sol.a = a1;
            sol.v = v1;

            fprintf('\n------------Life Cycle Problem Solved.------------\n');
        end
    end
end

