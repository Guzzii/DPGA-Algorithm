clear all
global n A b groups w ell_1 ell_2 delta nodes nGroups

nGroups = 10; % number of groups test with 10, 30, 50
delta = 1; % Huber loss function parameter
ell_1 = 1; % Coefficient for L-1 norm
ell_2 = 1; % Coefficient for group norm

result_dpga_100 = cell(5, 4);
k = 0;
for seed = 13:13
    for gsize = [100]
        for nodes = [100]
            for group_flag = [true] % true for same grouping 
                for net = [3] % 1 for tree, 2 for clique, 3 for smallworld
                    for edges = [0, 100, 300, 700, 4850]; % # of edges in addtion to small world
                        k = k + 1;
                        [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                            seed, group_flag, edges );
                        [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                            A, b, w, seed, group_flag, 0);
                        [x, relfun_dpga, infeas_dpga, walltime_dpga, normsqd, L] = ...
                            DPGA( nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
                            cvx_optval, xc, net, seed, group_flag, 1, 1);
                        
                        r = 3.5156*sqrt(4/size(E, 1)/min(diag(G)));
                        fun_bnd_cons = 4/r/min(eigs(G))*normsqd + 0.5*sum(L)*(xc)'*xc;
                        
                        result_dpga_100{k, 1} = relfun_dpga;
                        result_dpga_100{k, 2} = infeas_dpga;
                        result_dpga_100{k, 3} = walltime_dpga;
                        result_dpga_100{k, 4} = fun_bnd_cons/cvx_optval;
                    end
                end
            end
        end
    end
end

%%
result_dpga_500 = cell(5, 4);
k = 0;
for seed = 13:13
    for gsize = [100]
        for nodes = [500]
            for group_flag = [true] % true for same grouping 
                for net = [3] % 1 for tree, 2 for clique, 3 for smallworld
                    for edges = [500, 3500]; % # of edges in addtion to small world
                        k = k + 1;
                        [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                            seed, group_flag, edges );
                        [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                            A, b, w, seed, group_flag, 0);
                        [x, relfun_dpga, infeas_dpga, walltime_dpga, normsqd, L] = ...
                            DPGA( nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
                            cvx_optval, xc, net, seed, group_flag, 1, 1);
                        
                        r = 3.5156*sqrt(4/size(E, 1)/min(diag(G)));
                        fun_bnd_cons = 4/r/min(diag(G))*normsqd + 0.5*sum(L)*(xc)'*xc;
                        
                        result_dpga_500{k, 1} = relfun_dpga;
                        result_dpga_500{k, 2} = infeas_dpga;
                        result_dpga_500{k, 3} = walltime_dpga;
                        result_dpga_500{k, 4} = fun_bnd_cons/cvx_optval;
                    end
                end
            end
        end
    end
end
load result_dpga.mat
save result_dpga.mat result_dpga_10 result_dpga_50 result_dpga_100 result_dpga_500