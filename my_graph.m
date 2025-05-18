%% File Info.
%{
    my_graph.m
    ----------
    This code plots:
    (1) Value and policy functions from the simulation
    (2) Smoothed income profile Gₜ from data
    (3) Life-cycle consumption and wealth profiles (inlcuding multiple
    parameters and heat map)
%}

%% Graph class.
classdef my_graph
    methods(Static)

        %% Plot value and policy functions.
        function [] = plot_policy(par, sol, sim)
            astate = par.agrid;
            age = linspace(1, par.T, par.T);

            % Consumption policy
            figure;
            surf(age(1:5:end), astate(1:15:end), sol.c(1:15:end,1:5:end))
            xlabel('t','Interpreter','latex')
            ylabel('$a_t$','Interpreter','latex')
            zlabel('$c_t$','Interpreter','latex')
            title('Consumption Policy Function')

            % Saving policy
            figure;
            surf(age(1:5:end), astate(1:15:end), sol.a(1:15:end,1:5:end))
            xlabel('t','Interpreter','latex')
            ylabel('$a_t$','Interpreter','latex')
            zlabel('$a_{t+1}$','Interpreter','latex')
            title('Saving Policy Function')

            % Value function
            figure;
            surf(age(1:5:end), astate(1:15:end), sol.v(1:15:end,1:5:end))
            xlabel('t','Interpreter','latex')
            ylabel('$a_t$','Interpreter','latex')
            zlabel('$v_t(a_t,t)$','Interpreter','latex')
            title('Value Function')
        end

        %% Plot smoothed income profile Gt
        function [] = plot_gt(Gt_bin)
            Gt_smooth = movmean(Gt_bin.Gt, [1 1], 'Endpoints', 'shrink');

            figure;
            plot(Gt_bin.AgeBin, Gt_smooth, 'r-', 'LineWidth', 2);
            title("Gₜ: Smoothed Geometric Mean of Income by Age (All Observed)");
            xlabel("Age Bin (e.g., 25 = ages 25–29)");
            ylabel("Gₜ = exp(mean(log income))");
            grid on;
        end

        %% Plot life-cycle profiles of consumption and wealth
        function [] = plot_lifecycle_profiles(par, sim, model_name)
            T = par.T;
            age = (0:T-1)';

            % Compute average consumption and wealth over time
            avg_c = mean(sim.csim, 2, 'omitnan');
            avg_a = mean(sim.asim, 2, 'omitnan');

            % Plot average consumption profile
            figure;
            plot(age, avg_c, 'b-', 'LineWidth', 2);
            title(['Average Consumption by Age - ', model_name]);
            xlabel('Age');
            ylabel('Average Consumption');
            grid on;

            % Plot average wealth profile
            figure;
            plot(age, avg_a, 'g-', 'LineWidth', 2);
            title(['Average Wealth by Age - ', model_name]);
            xlabel('Age');
            ylabel('Average Wealth');
            grid on;
        end

        %% Average life-cycle profiles of consumption and wealth
        function [avg_c, avg_a] = average_profiles(par, sim)
            T = par.T;
            avg_c = nan(T,1);
            avg_a = nan(T,1);
            for t = 1:T
                avg_c(t) = mean(sim.csim(t, :), 'omitnan');
                avg_a(t) = mean(sim.asim(t, :), 'omitnan');
            end
        end

        %% Plot multiple life-cycle profiles for varying beta or gamma
        function [] = plot_multiple_profiles(par, profiles, param_values, param_name, var_type)
            figure;
            hold on;
            colors = lines(length(param_values));

            for i = 1:length(param_values)
                if ~isempty(profiles{i})
                    plot(0:par.T-1, profiles{i}, 'Color', colors(i,:), 'LineWidth', 2, ...
                        'DisplayName', sprintf('%s = %.2f', param_name, param_values(i)));
                end
            end

            hold off;
            legend('Location', 'best');
            title(['Life-Cycle ', var_type, ' Profiles (Varying ', param_name, ')']);
            xlabel('Age');
            ylabel(['Average ', var_type]);
            grid on;
        end

        %% Plot heatmap of average wealth
        function [] = plot_heatmap(avg_wealth_matrix, betas, gammas)
            figure;
            imagesc(gammas, betas, avg_wealth_matrix);
            colorbar;
            xlabel('$\gamma$ (Risk Aversion)', 'Interpreter', 'latex');
            ylabel('$\beta$ (Discount Factor)', 'Interpreter', 'latex');
            title('Average Wealth Heatmap', 'Interpreter', 'latex');
            set(gca, 'YDir', 'normal'); 
            xticks(gammas);
            yticks(betas);
        end
    end
end

