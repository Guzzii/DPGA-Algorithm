
% 100 nodes network
clear h
output = main({'dpga'}, 100);
result_dpga_100 = output.result_dpga;

f = @(x, y) y./x;

h(1) = figure;
semilogy(1:length(result_dpga_100{1, 1}), result_dpga_100{1, 1}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_100{2, 1}), result_dpga_100{2, 1}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{3, 1}), result_dpga_100{3, 1}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{4, 1}), result_dpga_100{4, 1}(1:end), '-r', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{5, 1}), result_dpga_100{5, 1}(1:end), '-g', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{1, 1}), f(1:length(result_dpga_100{1, 1}), result_dpga_100{1, 4}),...
    '-b', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('Circle', '200 Edges', '400 Edges', '800 Edges', 'Complete', 'Theoratical Bound');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Relative Suboptimality', 'FontSize', 18)
grid on
hold off

h(2) = figure;
semilogy(1:length(result_dpga_100{1, 2}), result_dpga_100{1, 2}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_100{2, 2}), result_dpga_100{2, 2}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{3, 2}), result_dpga_100{3, 2}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{4, 2}), result_dpga_100{4, 2}(1:end), '-r', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{5, 2}), result_dpga_100{5, 2}(1:end), '-g', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{1, 2}), f(1:length(result_dpga_100{1, 2}), 0.005*result_dpga_100{1, 4}),...
    '-b', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('Circle', '200 Edges', '400 Edges', '800 Edges', 'Complete', 'Theoratical Bound');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Infeasibility', 'FontSize', 18)
grid on
hold off

savefig(h, 'DPGA_100_Nodes.fig')

%%
% 100 nodes network
clear h
output = main({'dpga'}, 100);
result_dpga_100 = output.result_dpga;

f = @(x, y) y./x;

h(1) = figure;
semilogy(1:length(result_dpga_100{1, 1}), result_dpga_100{1, 1}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_100{2, 1}), result_dpga_100{2, 1}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{3, 1}), result_dpga_100{3, 1}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{4, 1}), result_dpga_100{4, 1}(1:end), '-r', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{5, 1}), result_dpga_100{5, 1}(1:end), '-g', 'LineWidth', 2);
xlim([0,30000])
plotLegend = legend('Circle', '200 Edges', '400 Edges', '800 Edges', 'Complete');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Relative Suboptimality', 'FontSize', 18)
grid on
hold off

h(2) = figure;
semilogy(1:length(result_dpga_100{1, 2}), result_dpga_100{1, 2}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_100{2, 2}), result_dpga_100{2, 2}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{3, 2}), result_dpga_100{3, 2}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{4, 2}), result_dpga_100{4, 2}(1:end), '-r', 'LineWidth', 2);
semilogy(1:length(result_dpga_100{5, 2}), result_dpga_100{5, 2}(1:end), '-g', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('Circle', '200 Edges', '400 Edges', '800 Edges', 'Complete');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Infeasibility', 'FontSize', 18)
grid on
hold off

savefig(h, 'DPGA_100_Nodes_no_Bnd.fig')

%%
clear h
output = main({'dpga'}, 100);
result_dpga_100 = output.result_dpga;
output = main({'dpga'}, 500);
result_dpga_500 = output.result_dpga;

h(1) = figure;
semilogy(1:length(result_dpga_100{2, 1}), result_dpga_100{2, 1}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_100{3, 1}), result_dpga_100{3, 1}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_500{1, 1}), result_dpga_500{1, 1}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_500{2, 1}), result_dpga_500{2, 1}(1:end), '-r', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('100 Nodes 200 Edges', '100 Nodes 400 Edges', '500 Nodes 1000 Edges',...
    '500 Nodes 4000 Edges');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Relative Suboptimality', 'FontSize', 18)
grid on
hold off

h(2) = figure;
semilogy(1:length(result_dpga_100{2, 2}), result_dpga_100{2, 2}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_100{3, 2}), result_dpga_100{3, 2}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_500{1, 2}), result_dpga_500{1, 2}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_500{2, 2}), result_dpga_500{2, 2}(1:end), '-r', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('100 Nodes 200 Edges', '100 Nodes 400 Edges', '500 Nodes 1000 Edges',...
    '500 Nodes 4000 Edges');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Infeasibility', 'FontSize', 18)
grid on
hold off

savefig(h, 'DPGA_Density_100_Nodes.fig')

%%
clear h
output = main({'dpga'}, 100);
result_dpga_100 = output.result_dpga;

f = @(x, y) y./x;

network = {'Circle', '200 Edges', '400 Edges', '800 Edges', 'Complete'};
ylab = {'Relative Suboptimality', 'Infeasibility'};

k=0;
for i = 1:2
    for j = 1:5
        k = k+1;
        h(k) = figure;
        semilogy(1:length(result_dpga_100{j, i}), result_dpga_100{j, i}(1:end), '-y', 'LineWidth', 2);
        hold on
        if i==1
            semilogy(1:length(result_dpga_100{j, i}), f(1:length(result_dpga_100{j, i}), ...
                result_dpga_100{1, 4}), '-m', 'LineWidth', 2);
        else
            semilogy(1:length(result_dpga_100{j, i}), f(1:length(result_dpga_100{j, i}), ...
                0.005*result_dpga_100{1, 4}), '-m', 'LineWidth', 2);
        end
        plotLegend = legend(network{j}, 'Theoratical Bound');
        set(plotLegend, 'FontSize', 18);
        xlabel('Iterations', 'FontSize', 18)
        ylabel(ylab{i}, 'FontSize', 18)
        grid on
        hold off
    end
end
savefig(h,'DPGA_100_Nodes_Bnd.fig')

%%
clear h
output = main({'dpga', 'pg_extra'}, 100);
result_dpga_100 = output.result_dpga;
result_pgextra_100 = output.result_pgextra;

xlab = {'Iterations', 'Communications', 'Walltime'};
ylab = {'Relative Suboptimality', 'Infeasibility'};

k=0;
for i = 1:2
    for j = 1:3
        k = k+1;
        h(k) = figure;
        if j==2
            semilogy([1:length(result_dpga_100{3, j})], result_dpga_100{3, i}(1:end), '-y', 'LineWidth', 2);
            hold on
            semilogy(2*[1:length(result_pgextra_100{1, j})], result_pgextra_100{1, i}(1:end),...
                '-m', 'LineWidth', 2);
        else
            semilogy(1:length(result_dpga_100{3, j}), result_dpga_100{3, i}(1:end), '-y', 'LineWidth', 2);
            hold on
            semilogy(1:length(result_pgextra_100{1, j}), result_pgextra_100{1, i}(1:end),...
                '-m', 'LineWidth', 2);
        end
        xlim([0, 30000])
        plotLegend = legend('DPGA', 'PG-EXTRA');
        set(plotLegend, 'FontSize', 18);
        xlabel(xlab{j}, 'FontSize', 18)
        ylabel(ylab{i}, 'FontSize', 18)
        grid on
        hold off
    end
end
savefig(h,'DPGA_vs_PG_EXTRA_100_Nodes.fig')

%%
clear h
output = main({'dpga', 'sdpga'}, 100);
result_dpga_100 = output.result_dpga;
result_sdpga_100 = output.result_sdpga;

step = 1000;
xlab = {'Iterations'};
ylab = {'Relative Suboptimality', 'Infeasibility'};

for i = 1:2
    h(i) = figure;
    semilogy(1:length(result_dpga_100{2, i}), result_dpga_100{2, i}(1:end), '-o', 'LineWidth', 2);
    hold on
    semilogy(1:length(result_sdpga_100{1, i}), result_sdpga_100{1, i}(1:end),...
        '-o', 'LineWidth', 2);
    semilogy(1:length(result_sdpga_100{2, i}), result_sdpga_100{2, i}(1:end),...
        '-o', 'LineWidth', 2);
    semilogy(1:length(result_sdpga_100{3, i}), result_sdpga_100{3, i}(1:end),...
        '-o', 'LineWidth', 2);
    xlim([0, 30000])
    plotLegend = legend('DPGA', 'SDPGA with \sigma = 1', 'SDPGA with \sigma = 0.1',...
        'SDPGA with \sigma = 0.01');
    set(plotLegend, 'FontSize', 18);
    xlabel(xlab{1}, 'FontSize', 18)
    ylabel(ylab{i}, 'FontSize', 18)
    grid on
    hold off
end
savefig(h,'DPGA_vs_SDPGA_100.fig')
