function [ output ] = main(algo_list, nodes)
%The main wrapper function for producing plots

    global n A b groups w ell_1 ell_2 delta nGroups

    nGroups = 10; % number of groups test with 10, 30, 50
    delta = 1; % Huber loss function parameter
    ell_1 = 1; % Coefficient for L-1 norm
    ell_2 = 1; % Coefficient for group norm
    seed = 13;
    gsize = 100;
    group_flag = true;
    net = 3; % 1 for tree, 2 for clique, 3 for smallworld
    
    output = [];
    
    if nodes == 10
        edge_list = [0, 10, 20, 35];
        edge_list_pg = 10;
        edge_list_sdpga = 10;
    elseif nodes == 100
        edge_list = [0, 100, 300, 700, 4850];
        edge_list_pg = 100;
        edge_list_sdpga = 100;
    elseif nodes == 50
        edge_list = [50, 100];
    elseif nodes == 500
        edge_list = [500, 3500];
    end
    
    if any(ismember(algo_list, 'dpga'))
        
        result_dpga = cell(5, 4);
        k = 0;
        for edges = edge_list; % # of edges in addtion to small world
            k = k + 1;
            [n, E, G, A, b, groups, w] = Data_Generate(nodes, gsize,...
                nGroups, net, seed, group_flag, edges);
            [cvx_optval, xc] = CVX_Central(nodes, gsize, n, nGroups, groups,...
                E, G, A, b, w, seed, group_flag, 0);
            [~, relfun_dpga, infeas_dpga, walltime_dpga, normsqd, L] = ...
                DPGA(nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
                cvx_optval, xc, net, seed, group_flag, 1800, 1);

            r = 3.5156*sqrt(4/size(E, 1)/min(diag(G)));
            fun_bnd_cons = 4/r/min(eigs(G))*normsqd + 0.5*sum(L)*(xc)'*xc;

            result_dpga{k, 1} = relfun_dpga;
            result_dpga{k, 2} = infeas_dpga;
            result_dpga{k, 3} = walltime_dpga;
            result_dpga{k, 4} = fun_bnd_cons/cvx_optval;
        end
        output.result_dpga = result_dpga;
        
    end
        
    if any(ismember(algo_list, 'pg_extra'))
        
        result_pgextra = cell(5, 3);
        k = 0;
        for edges = edge_list_pg; % # of edges in addtion to small world
            k = k + 1;
            [n, E, G, A, b, groups, w] = Data_Generate(nodes, gsize,...
                nGroups, net, seed, group_flag, edges);
            [cvx_optval, xc] = CVX_Central(nodes, gsize, n, nGroups, groups,...
                E, G, A, b, w, seed, group_flag, 0);
            [~, relfun_pgextra, infeas_pgextra, walltime_pgextra] = ...
                PG_EXTRA( nodes, gsize, n, nGroups, groups, E, G, A, b, w, ...
                cvx_optval, net, seed, group_flag);

            result_pgextra{k, 1} = relfun_pgextra;
            result_pgextra{k, 2} = infeas_pgextra;
            result_pgextra{k, 3} = walltime_pgextra;
        end
        output.result_pgextra = result_pgextra;

    end
    
    if any(ismember(algo_list, 'sdpga'))
        
        result_sdpga = cell(5, 3);
        k = 0;
        for edges = edge_list_sdpga; % # of edges in addtion to small world
            for sigma = [100, 10, 1] % sdv for the noise
                k = k + 1;
                [n, E, G, A, b, groups, w] = Data_Generate(nodes, gsize,...
                    nGroups, net, seed, group_flag, edges);
                [cvx_optval, xc] = CVX_Central(nodes, gsize, n, nGroups,...
                    groups, E, G, A, b, w, seed, group_flag, 0);
                [~, relfun_sdpga, infeas_sdpga, walltime_sdpga] = ...
                    SDPGA(nodes, gsize, n, nGroups, groups, E, G, A, b, w,...
                    cvx_optval, net, seed, group_flag, 600, sigma);
                
                result_sdpga{k, 1} = relfun_sdpga;
                result_sdpga{k, 2} = infeas_sdpga;
                result_sdpga{k, 3} = walltime_sdpga;
            end
        end
        output.result_sdpga = result_sdpga;
                        
    end
    
end
