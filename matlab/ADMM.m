function[x, Relfun, Infeas, Walltime] = ADMM(nodes, gsize, n, nGroups, groups, E, G, A, b, w, cvx_optval, net, seed, group_flag)

global ell_1 ell_2 delta c

if net==1
    if (gsize==100&&(nodes==5||nodes==10))
        c = 4; % 5&10, 10, 100, tree
    else
        c = 0.54;
    end
else
    if (gsize==100&&(nodes==5||nodes==10))
        c = 0.75; % 5&10, 10, 100, clique
    else
        c = 0.42;
    end
end

arc=size( E, 1); % Number of arcs in the graph
x = cell(nodes,1);
p = cell(nodes, 1);
for i=1:nodes
    x{i} = zeros(n, 1);
    p{i} = zeros(n, 1);
end
%ADMM SOLUTION
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
            y = G_x(G, x, 2);
        end
        % Solve subproblems
        for i=1:nodes
            prox_count = prox_count + 1;
            x{i} = prox1( x{i}, p, y, G, i, 1 );
        end
        % Update dual variables
        y = G_x(G, x, 2);
        for i=1:nodes
            p{i} = p{i}+c*y{i};
        end
        % Compute the objective function value at current iterates
        val1 = 0;
        for i=1:nodes
            for k=1:nGroups
                val1 = val1+ell_2/nodes*w(k)*norm(x{i}(groups{k}));
            end
            val1 = val1+ell_1/nodes*norm(x{i},1) + sum( huberloss( A{i}*x{i}-b{i}, delta ) );
        end
        %Compute infeasiblity
        dist=0;
        for i=1:arc
            arcList=find(E(i,:)~=0);
            dist=max( dist, norm(x{arcList(1)}-x{arcList(2)}, 2)/sqrt(n) );
        end
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
            y = G_x(G, x, 2);
        end
        % Solve subproblems
        for i=1:nodes
            prox_count = prox_count + 1;
            x{i} = prox1( x{i}, p, y, G, i, 2 );
        end
        % Update dual variables
        y = G_x(G, x, 2);
        for i=1:nodes
            p{i} = p{i}+c*y{i};
        end
        % Compute the objective function value at current iterates
        val1 = 0;
        for i=1:nodes
            for k=1:nGroups
                val1 = val1+ell_2/nodes*w(k)*norm(x{i}(groups{i, k}));
            end
            val1 = val1+ell_1/nodes*norm(x{i},1) + sum( huberloss( A{i}*x{i}-b{i}, delta ) );
        end
        %Compute infeasiblity
        dist=0;
        for i=1:arc
            arcList=find(E(i,:)~=0);
            dist=max( dist, norm(x{arcList(1)}-x{arcList(2)}, 2)/sqrt(n) );
        end
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
t = toc

if group_flag
    if net == 1
        row_name = ['ADMM_tree','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    else
        row_name = ['ADMM_clique','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    end
else
    if net == 1
        row_name = ['ADMM_tree','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    else
        row_name = ['ADMM_clique','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    end
end
output_admm = dataset({[relfun, dist, t, outer_iter, prox_count/nodes], 'Rel', 'ConsensusViolation',...
    'WallTime', 'Iterations', 'ProxOperations'}, 'ObsNames', {row_name});
if isempty(dir('Result.csv'))
    output = output_admm;
    export(output, 'File', 'Result.csv', 'Delimiter', ',')
else
    output_old = dataset('File', 'Result.csv', 'Delimiter', ',', 'ReadObsNames', true);
    if ismember(row_name, output_old.Properties.ObsNames)
        output_old(row_name, :) = output_admm;
        output_new = output_old;
    else
        output_new = [output_old; output_admm];
    end
    export(output_new, 'File', 'Result.csv', 'Delimiter', ',')
end

end

