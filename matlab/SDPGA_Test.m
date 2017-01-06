clear all
global n A b groups w ell_1 ell_2 delta nodes nGroups

nGroups = 10; % number of groups test with 10, 30, 50
delta = 1; % Huber loss function parameter
ell_1 = 1; % Coefficient for L-1 norm
ell_2 = 1; % Coefficient for group norm

result_sdpga_10 = cell(5, 3);
k = 0;
for seed = 13:13
    for gsize = [100]
        for nodes = [10]
            for group_flag = [true] % true for same grouping 
                for net = [3] % 1 for tree, 2 for clique, 3 for smallworld
                    for edges = 10 % # of edges in addtion to small world
                        for sigma = [1, 0.1, 0.01] % sdv for the noise
                            k = k + 1;
                            [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                                seed, group_flag, edges );
                            [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                                A, b, w, seed, group_flag, 0);

                            [xsii, relfun_sdpga_ii, infeas_sdpga_ii, walltime_sdpga_ii, normsqds, Ls] = ...
                                SDPGA_II( nodes, gsize, n, nGroups, groups, ...
                                E, G, A, b, w, cvx_optval, net, seed, group_flag, 600, sigma);

                            r = 25.5156*sqrt(4/size(E, 1)/min(diag(G)));
                            fun_bnd_cons = 4/r/min(diag(G))*normsqds + 0.5*sum(Ls)*(xc)'*xc;

                            result_sdpga_10{k, 1} = relfun_sdpga_ii;
                            result_sdpga_10{k, 2} = infeas_sdpga_ii;
                            result_sdpga_10{k, 3} = walltime_sdpga_ii;
                            result_sdpga_10{k, 4} = fun_bnd_cons/cvx_optval;
                        end
                    end
                end
            end
        end
    end
end

result_sdpga_100 = cell(5, 3);
k = 0;
for seed = 13:13
    for gsize = [100]
        for nodes = [100]
            for group_flag = [true] % true for same grouping 
                for net = [3] % 1 for tree, 2 for clique, 3 for smallworld
                    for edges = 100 % # of edges in addtion to small world
                        for sigma = [1, 0.1, 0.01] % sdv for the noise
                            k = k + 1;
                            [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                                seed, group_flag, edges );
                            [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                                A, b, w, seed, group_flag, 0);

                            [xsii, relfun_sdpga_ii, infeas_sdpga_ii, walltime_sdpga_ii, normsqds, Ls] = ...
                                SDPGA_II( nodes, gsize, n, nGroups, groups, ...
                                E, G, A, b, w, cvx_optval, net, seed, group_flag, 600, sigma);

                            r = 25.5156*sqrt(4/size(E, 1)/min(diag(G)));
                            fun_bnd_cons = 4/r/min(diag(G))*normsqds + 0.5*sum(Ls)*(xc)'*xc;

                            result_sdpga_100{k, 1} = relfun_sdpga_ii;
                            result_sdpga_100{k, 2} = infeas_sdpga_ii;
                            result_sdpga_100{k, 3} = walltime_sdpga_ii;
                            result_sdpga_100{k, 4} = fun_bnd_cons/cvx_optval;
                        end
                    end
                end
            end
        end
    end
end

save result_sdpga.mat result_sdpga_10 result_sdpga_100