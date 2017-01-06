function[x, Relfun, Infeas, Walltime] = SADMM(nodes, gsize, n, nGroups, groups, E, G, A, b, w, cvx_optval, net, seed, group_flag)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global ell_1 ell_2 delta c

if net==1
    if gsize==100
        c = 1; % 5&10, 10, 100, tree
    else
        c = 0.36;
    end
else
    if gsize==100
        c = 0.5; % 5&10, 10, 100, clique
    else
        c = 0.28;
    end
end


arc=size( E, 1); % Number of arcs in the graph
x = cell(nodes,1);
p = cell(nodes, 1);
for i=1:nodes
    x{i} = zeros(n, 1);
    p{i} = zeros(n, 1);
end
y = x;
r = x;
pbar = p;
%SADMM SOLUTION
tic

tol1 = 1e-3;
tol2 = 1e-4;
opt_flag = 1;
prox_count = 0;
outer_iter = 0;
Relfun = [];
Infeas = [];
Walltime = [];

if group_flag
    while(opt_flag)
        if outer_iter==0
            s = G_x(G, x, 2);
            sbar = G_x(G, y, 2);
        end
        % Solve subproblems
        for i=1:nodes
            prox_count = prox_count + 2;
            x{i} = prox2( x{i}, y{i}, r{i}, p, s, G, i, 3 );
            y{i} = prox2( x{i}, y{i}, r{i}, pbar, sbar, G, i, 4 );
        end
        % Update dual variables
        s = G_x(G, x, 2);
        sbar = G_x(G, y, 2);
        for i=1:nodes
            p{i} = p{i}+c*s{i};
            pbar{i} = pbar{i}+c*sbar{i};
            r{i} = r{i} + (x{i}-y{i})/2;
        end
        % Compute the objective function value at current iterates
        val1 = 0;
        for i=1:nodes
            for k=1:nGroups
                val1 = val1+ell_2/nodes*w(k)*norm((x{i}(groups{k})+y{i}(groups{k}))/2);
            end
            val1 = val1+ell_1/nodes*norm((x{i}+y{i})/2,1) + sum( huberloss( A{i}*(x{i}+y{i})/2-b{i}, delta ) );
        end
        %Compute infeasiblity
        dist1 = 0;
        dist2 = 0;
        dist3 = 0;
        for i=1:arc
            arcList = find(E(i,:)~=0);
            dist1 = max( dist1, norm(x{arcList(1)}-x{arcList(2)}, 2)/sqrt(n) );
            dist3 = max( dist3, norm(y{arcList(1)}-y{arcList(2)}, 2)/sqrt(n) );
        end
        for i=1:nodes
            dist2 = max( dist2, norm(x{i} -y{i})/sqrt(n) );
        end
        dist = max( [dist1, dist2, dist3] );
        relfun = abs( cvx_optval - val1)/cvx_optval;

        %Check stopping criteria
        if relfun < tol1&&dist<tol2
            opt_flag = 0;
        end
        
        Relfun = [Relfun, relfun];
        Infeas = [Infeas, dist];

        t_break = toc;
        Walltime = [Walltime, t_break];
        
        if t_break>1800
            opt_flag = 0;
        end

        outer_iter = outer_iter+1;
        disp( ['Out_Ite: ', num2str(outer_iter), ' -- Obj: ', num2str( val1 ), ' -- |xi=xj|: ', num2str( dist )] );
    end
else
    while(opt_flag)
        if outer_iter==0
            s = G_x(G, x, 2);
            sbar = G_x(G, y, 2);
        end
        % Solve subproblems
        for i=1:nodes
            prox_count = prox_count + 2;
            x{i} = prox2( x{i}, y{i}, r{i}, p, s, G, i, 5 );
            y{i} = prox2( x{i}, y{i}, r{i}, pbar, sbar, G, i, 4 );
        end
        % Update dual variables
        s = G_x(G, x, 2);
        sbar = G_x(G, y, 2);
        for i=1:nodes
            p{i} = p{i}+c*s{i};
            pbar{i} = pbar{i}+c*sbar{i};
            r{i} = r{i} + (x{i}-y{i})/2;
        end
        % Compute the objective function value at current iterates
        val1 = 0;
        for i=1:nodes
            for k=1:nGroups
                val1 = val1+ell_2/nodes*w(k)*norm((x{i}(groups{i, k})+y{i}(groups{i, k}))/2);
            end
            val1 = val1+ell_1/nodes*norm((x{i}+y{i})/2,1) + sum( huberloss( A{i}*(x{i}+y{i})/2-b{i}, delta ) );
        end
        %Compute infeasiblity
        dist1 = 0;
        dist2 = 0;
        dist3 = 0;
        for i=1:arc
            arcList = find(E(i,:)~=0);
            dist1 = max( dist1, norm(x{arcList(1)}-x{arcList(2)}, 2)/sqrt(n) );
            dist3 = max( dist3, norm(y{arcList(1)}-y{arcList(2)}, 2)/sqrt(n) );
        end
        for i=1:nodes
            dist2 = max( dist2, norm(x{i} -y{i})/sqrt(n) );
        end
        dist = max( [dist1, dist2, dist3] );
        relfun = abs( cvx_optval - val1)/cvx_optval;

        %Check stopping criteria
        if relfun < tol1&&dist<tol2
            opt_flag = 0;
        end
        
        Relfun = [Relfun, relfun];
        Infeas = [Infeas, dist];

        t_break = toc;
        Walltime = [Walltime, t_break];
        
        if t_break>1800
            opt_flag = 0;
        end

        outer_iter = outer_iter+1;
        disp( ['Out_Ite: ', num2str(outer_iter), ' -- Obj: ', num2str( val1 ), ' -- |xi=xj|: ', num2str( dist )] );
    end
end
t=toc

if group_flag
    if net == 1
        row_name = ['SADMM_tree','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    else
        row_name = ['SADMM_clique','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    end
else
    if net == 1
        row_name = ['SADMM_tree','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    else
        row_name = ['SADMM_clique','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    end
end
output_sadmm = dataset({[relfun, dist, t, outer_iter, prox_count/nodes], 'Rel', 'ConsensusViolation',...
    'WallTime', 'Iterations', 'ProxOperations'}, 'ObsNames', {row_name});
if isempty(dir('Result.csv'))
    output = output_sadmm;
    export(output, 'File', 'Result.csv', 'Delimiter', ',')
else
    output_old = dataset('File', 'Result.csv', 'Delimiter', ',', 'ReadObsNames', true);
    if ismember(row_name, output_old.Properties.ObsNames)
        output_old(row_name, :) = output_sadmm;
        output_new = output_old;
    else
        output_new = [output_old; output_sadmm];
    end
    export(output_new, 'File', 'Result.csv', 'Delimiter', ',')
end
end

