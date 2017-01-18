% Script generates plots, each section produces one set of figures

%% DPGA plots with 10 nodes 
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

%% DPGA density plots: 10 vs 50 nodes
%
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

%% DPGA vs SDPGA with 10 nodes
%
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
