% Script generates plots when nodes = 10 and 50
clear h
output = main({'dpga'}, 10);
result_dpga_10 = output.result_dpga;

f = @(x, y) y./x;

h(1) = figure;
semilogy(1:length(result_dpga_10{1, 1}), result_dpga_10{1, 1}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_10{2, 1}), result_dpga_10{2, 1}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{3, 1}), result_dpga_10{3, 1}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{4, 1}), result_dpga_10{4, 1}(1:end), '-r', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{1, 1}), f(1:length(result_dpga_10{1, 1}), result_dpga_10{1, 4}),...
    '-b', 'LineWidth', 2);
xlim([0, 10000])
plotLegend = legend('Circle', '20 Edges', '30 Edges', 'Complete', 'Theoratical Bound');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Relative Suboptimality', 'FontSize', 18)
grid on
hold off

h(2) = figure;
semilogy(1:length(result_dpga_10{1, 2}), result_dpga_10{1, 2}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_10{2, 2}), result_dpga_10{2, 2}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{3, 2}), result_dpga_10{3, 2}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{4, 2}), result_dpga_10{4, 2}(1:end), '-r', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{1, 2}), f(1:length(result_dpga_10{1, 2}), 0.5*result_dpga_10{1, 4}),...
    '-b', 'LineWidth', 2);
xlim([0, 10000])
plotLegend = legend('Circle', '20 Edges', '30 Edges', 'Complete', 'Theoratical Bound');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Infeasibility', 'FontSize', 18)
grid on
hold off

savefig(h, 'DPGA_10_Nodes.fig')

%%
% 100 nodes network
clear h
output = main({'dpga'}, 10);
result_dpga_10 = output.result_dpga;

f = @(x, y) y./x;

h(1) = figure;
semilogy(1:length(result_dpga_10{1, 1}), result_dpga_10{1, 1}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_10{2, 1}), result_dpga_10{2, 1}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{3, 1}), result_dpga_10{3, 1}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{4, 1}), result_dpga_10{4, 1}(1:end), '-r', 'LineWidth', 2);
xlim([0, 10000])
plotLegend = legend('Circle', '20 Edges', '30 Edges', 'Complete');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Relative Suboptimality', 'FontSize', 18)
grid on
hold off

h(2) = figure;
semilogy(1:length(result_dpga_10{1, 2}), result_dpga_10{1, 2}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_10{2, 2}), result_dpga_10{2, 2}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{3, 2}), result_dpga_10{3, 2}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_10{4, 2}), result_dpga_10{4, 2}(1:end), '-r', 'LineWidth', 2);
xlim([0, 10000])
plotLegend = legend('Circle', '20 Edges', '30 Edges', 'Complete');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Infeasibility', 'FontSize', 18)
grid on
hold off

savefig(h, 'DPGA_10_Nodes_no_Bnd.fig')

%%
clear h
output = main({'dpga'}, 10);
result_dpga_10 = output.result_dpga;
output = main({'dpga'}, 50);
result_dpga_50 = output.result_dpga;

h(1) = figure;
semilogy(1:length(result_dpga_10{2, 1}), result_dpga_10{2, 1}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_10{3, 1}), result_dpga_10{3, 1}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_50{1, 1}), result_dpga_50{1, 1}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_50{2, 1}), result_dpga_50{2, 1}(1:end), '-r', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('10 Nodes 20 Edges', '10 Nodes 30 Edges', '50 Nodes 100 Edges',...
    '50 Nodes 150 Edges');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Relative Suboptimality', 'FontSize', 18)
grid on
hold off

h(2) = figure;
semilogy(1:length(result_dpga_10{2, 2}), result_dpga_10{2, 2}(1:end), '-y', 'LineWidth', 2);
hold on
semilogy(1:length(result_dpga_10{3, 2}), result_dpga_10{3, 2}(1:end),  '-m', 'LineWidth', 2);
semilogy(1:length(result_dpga_50{1, 2}), result_dpga_50{1, 2}(1:end), '-c', 'LineWidth', 2);
semilogy(1:length(result_dpga_50{2, 2}), result_dpga_50{2, 2}(1:end), '-r', 'LineWidth', 2);
xlim([0, 30000])
plotLegend = legend('10 Nodes 20 Edges', '10 Nodes 30 Edges', '50 Nodes 100 Edges',...
    '50 Nodes 150 Edges');
set(plotLegend, 'FontSize', 18);
xlabel('Iterations', 'FontSize', 18)
ylabel('Infeasibility', 'FontSize', 18)
grid on
hold off

savefig(h, 'DPGA_Density_10_Nodes.fig')

%%
clear h

output = main({'dpga'}, 10);
result_dpga_10 = output.result_dpga;

f = @(x, y) y./x;

network = {'Circle', '20 Edges', '30 Edges', 'Complete'};
ylab = {'Relative Suboptimality', 'Infeasibility'};

k=0;
for i = 1:2
    for j = 1:4
        k = k+1;
        h(k) = figure;
        semilogy(1:length(result_dpga_10{j, i}), result_dpga_10{j, i}(1:end), '-y', 'LineWidth', 2);
        hold on
        semilogy(1:length(result_dpga_10{j, i}), f(1:length(result_dpga_10{j, i}), result_dpga_10{1, 4}),...
        '-m', 'LineWidth', 2);
        xlim([0, 10000])
        plotLegend = legend(network{j}, 'Theoratical Bound');
        set(plotLegend, 'FontSize', 18);
        xlabel('Iterations', 'FontSize', 18)
        ylabel(ylab{i}, 'FontSize', 18)
        grid on
        hold off
    end
end
savefig(h,'DPGA_10_Nodes_Bnd.fig')

%%
clear h
output = main({'dpga', 'pg_extra'}, 10);
result_dpga_10 = output.result_dpga;
result_pgextra_10 = output.result_pgextra;


step = 1000;
xlab = {'Iterations', 'Communications', 'Walltime'};
ylab = {'Relative Suboptimality', 'Infeasibility'};

k=0;
for i = 1:2
    for j = 1:3
        k = k+1;
        h(k) = figure;
        if j==2
            semilogy([1:length(result_dpga_10{3, j})], result_dpga_10{3, i}(1:end), '-y', 'LineWidth', 2);
            hold on
            semilogy(2*[1:length(result_pgextra_10{1, j})], result_pgextra_10{1, i}(1:end),...
                '-m', 'LineWidth', 2);
        else
            semilogy(1:length(result_dpga_10{3, j}), result_dpga_10{3, i}(1:end), '-y', 'LineWidth', 2);
            hold on
            semilogy(1:length(result_pgextra_10{1, j}), result_pgextra_10{1, i}(1:end),...
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
savefig(h,'DPGA_vs_PG_EXTRA_10_Nodes.fig')

%%
clear h
output = main({'dpga', 'sdpga'}, 10);
result_dpga_10 = output.result_dpga;
result_sdpga_10 = output.result_sdpga;

step = 1000;
xlab = {'Iterations'};
ylab = {'Relative Suboptimality', 'Infeasibility'};

for i = 1:2
    h(i) = figure;
    semilogy(1:length(result_dpga_10{2, i}), result_dpga_10{2, i}(1:end), '-y', 'LineWidth', 2);
    hold on
    semilogy(1:length(result_sdpga_10{1, i}), result_sdpga_10{1, i}(1:end),...
        '-m', 'LineWidth', 2);
    semilogy(1:length(result_sdpga_10{2, i}), result_sdpga_10{2, i}(1:end),...
        '-c', 'LineWidth', 2);
    semilogy(1:length(result_sdpga_10{3, i}), result_sdpga_10{3, i}(1:end),...
        '-r', 'LineWidth', 2);
    plotLegend = legend('DPGA', 'SDPGA with \sigma = 100', 'SDPGA with \sigma = 10',...
        'SDPGA with \sigma = 1');
    set(plotLegend, 'FontSize', 18);
    xlabel(xlab{1}, 'FontSize', 18)
    ylabel(ylab{i}, 'FontSize', 18)
    grid on
    hold off
end
savefig(h,'DPGA_vs_SDPGA_10.fig')
