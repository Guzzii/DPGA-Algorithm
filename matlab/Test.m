clear all
global n A b groups w ell_1 ell_2 delta nodes nGroups

nGroups = 10; % number of groups test with 10, 30, 50
delta = 1; % Huber loss function parameter
ell_1 = 1; % Coefficient for L-1 norm
ell_2 = 1; % Coefficient for group norm

result_dpga_100 = cell(5, 4);
k = 1;
for seed = 13:13
    for gsize = [100]
        for nodes = [100]
            for group_flag = [true] % true for same grouping 
                for net = [3] % 1 for tree, 2 for clique, 3 for smallworld
                    for edges = [100]; % # of edges in addtion to small world
                        k = k + 1;
                        [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                            seed, group_flag, edges );
                        [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                            A, b, w, seed, group_flag, 0);
                        [xii, relfun_dpga_ii, infeas_dpga_ii, walltime_dpga_ii, normsqd, L] = ...
                            DPGA_II( nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
                            cvx_optval, net, seed, group_flag);
                        
                        r = 2.5156*sqrt(4/size(E, 1)/min(diag(G)));
                        fun_bnd_cons = 4/r/min(diag(G))*normsqd + 0.5*sum(1/L)*(xc)'*xc;
                        
                        result_dpga_100{k, 1} = relfun_dpga_ii;
                        result_dpga_100{k, 2} = infeas_dpga_ii;
                        result_dpga_100{k, 3} = walltime_dpga_ii;
                        result_dpga_100{k, 4} = fun_bnd_cons;
                    end
                end
            end
        end
    end
end

%%
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
                    for edges = 20 % # of edges in addtion to small world
                        for sigma = [1, 0.1, 0.01] % sdv for the noise
                            k = k + 1;
                            [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                                seed, group_flag, edges );
                            [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                                A, b, w, seed, group_flag, 0);

                            [xsii, relfun_sdpga_ii, infeas_sdpga_ii, walltime_sdpga_ii, normsqds, Ls] = ...
                                SDPGA_II( nodes, gsize, n, nGroups, groups, ...
                                E, G, A, b, w, cvx_optval, xc, net, seed, group_flag, 600, sigma);

                            r = 25.5156*sqrt(4/size(E, 1)/min(diag(G)));
                            fun_bnd_cons = 4/r/min(diag(G))*normsqds + 0.5*sum(Ls)*(xc)'*xc;

                            result_sdpga_10{k, 1} = relfun_sdpga_ii;
                            result_sdpga_10{k, 2} = infeas_sdpga_ii;
                            result_sdpga_10{k, 3} = walltime_sdpga_ii;
                            result_sdpga_10{k, 4} = fun_bnd_cons;
                        end
                    end
                end
            end
        end
    end
end

save result_sdpga.mat result_sdpga_10

%%
clear all
global n A b groups w ell_1 ell_2 delta nodes nGroups

nGroups = 10; % number of groups test with 10, 30, 50
delta = 1; % Huber loss function parameter
ell_1 = 1; % Coefficient for L-1 norm
ell_2 = 1; % Coefficient for group norm

result_pgextra_10 = cell(5, 3);
k = 0;
for seed = 13:13
    for gsize = [100]
        for nodes = [10]
            for group_flag = [true] % true for same grouping 
                for net = [3] % 1 for tree, 2 for clique, 3 for smallworld
                    for edges = 20; %[0, 100, 300, 700, 4850]; % # of edges in addtion to small world
                        k = k + 1;
                        [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
                            seed, group_flag, edges );
                        [ cvx_optval, xc ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
                            A, b, w, seed, group_flag, 0);
                        [xii, relfun_pgextra, infeas_pgextra, walltime_pgextra] = ...
                            PG_EXTRA( nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
                            cvx_optval, net, seed, group_flag);
                                                
                        result_pgextra_10{k, 1} = relfun_pgextra;
                        result_pgextra_10{k, 2} = infeas_pgextra;
                        result_pgextra_10{k, 3} = walltime_pgextra;
                    end
                end
            end
        end
    end
end
save result_pg_extra.mat result_pgextra_10
%%
% %%
% for seed = 13:13
%     for gsize = [100]
%         for nodes = [5]
%             for group_flag = [true] % true for same grouping 
%                 for net = [2] % 1 for tree, 2 for clique, 3 for smallworld
%                     for edges = 20:20; % # of edges in addtion to the cycle
%                     [ n, E, G, A, b, groups, w ] = Data_Generate(nodes, gsize, nGroups, net,...
%                         seed, group_flag, edges );
%                     [ cvx_optval ] = CVX_Central(nodes, gsize, n, nGroups, groups, E, G, ...
%                         A, b, w, seed, group_flag, 0);
% %                     [xi, relfun_dpga_i, infeas_dpga_i, walltime_dpga_i] = DPGA_I(nodes, gsize, n, nGroups, ...
% %                         groups, E, G, A, b, w, cvx_optval, net, seed, group_flag);
%                     [xii, relfun_dpga_ii, infeas_dpga_ii, walltime_dpga_ii] = DPGA_II( nodes, gsize, n, ...
%                         nGroups, groups, E, G, A, b, w, cvx_optval, net, seed, group_flag);
% %                     [ xpg, relfun_pg, infeas_pg, walltime_pg] = PG_EXTRA(nodes, gsize, n, nGroups, groups,...
% %                         E, G, A, b, w, cvx_optval, net, seed, group_flag);
% %                     [ xad, relfun_ad, infeas_ad, walltime_ad] = ADMM(nodes, gsize, n, nGroups, groups, E, G,...
% %                         A, b, w, cvx_optval, net, seed, group_flag);
% %                     [ xsad, relfun_sad, infeas_sad, walltime_sad] = SADMM(nodes, gsize, n, nGroups, groups,...
% %                         E, G, A, b, w, cvx_optval, net, seed, group_flag);
%                     filename = ['plotData_', num2str(seed), '_', num2str(gsize), '_', num2str(nodes), '_',...
%                         num2str(group_flag), '_', num2str(net)];
%                     save(filename, 'relfun_dpga_i', 'relfun_dpga_ii', 'relfun_pg', 'relfun_ad', 'relfun_sad', ...
%                        'infeas_dpga_i', 'infeas_dpga_ii','infeas_pg', 'infeas_ad','infeas_sad',...
%                        'walltime_dpga_i','walltime_dpga_ii','walltime_pg','walltime_ad','walltime_sad')
%                     end
%                 end
%             end
%         end
%     end
% end