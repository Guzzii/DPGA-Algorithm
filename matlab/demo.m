clear all
global n A b groups w ell_1 ell_2 delta nodes nGroups


delta = 1; % Huber loss function parameter
ell_1 = 1; % Coefficient for L-1 norm
ell_2 = 1; % Coefficient for group norm
nGroups = 10; % number of groups
gsize = 100; % group size for each group

seed = 13;
nodes = 10;
group_flag = true; % true for same grouping
net = 3; % 1 for tree, 2 for clique, 3 for smallworld
edges = 20; % # of edges in addtion to the cycle in small world

[n, E, G, A, b, groups, w] = Data_Generate(nodes, gsize,...
    nGroups, net, seed, group_flag, edges);
[cvx_optval, xc] = CVX_Central(nodes, gsize, n, nGroups, groups,...
    E, G, A, b, w, seed, group_flag, 0);
% DPGA
[xdp, relfun_dpga, infeas_dpga, walltime_dpga, normsqd, L] = ...
    DPGA(nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
    cvx_optval, xc, net, seed, group_flag, 1800, 1);
% PG_EXTRA
[xpg, relfun_pgextra, infeas_pgextra, walltime_pgextra] = ...
    PG_EXTRA( nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
    cvx_optval, net, seed, group_flag);
% SDPGA
[xsdp, relfun_sdpga, infeas_sdpga, walltime_sdpga] = ...
    SDPGA(nodes, gsize, n, nGroups, groups, E, G, A, b, w,...
    cvx_optval, net, seed, group_flag, 600, 10);
%vADMM
[xad, relfun_ad, infeas_ad, walltime_ad] = ...
    ADMM(nodes, gsize, n, nGroups, groups, E, G,...
    A, b, w, cvx_optval, net, seed, group_flag);
% SADMM
[xsad, relfun_sad, infeas_sad, walltime_sad] = ...
    SADMM(nodes, gsize, n, nGroups, groups,...
    E, G, A, b, w, cvx_optval, net, seed, group_flag);

