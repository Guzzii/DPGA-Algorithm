function [x, Relfun, Infeas, Walltime] = PG_EXTRA( nodes, gsize, n, nGroups, groups, E, G, A, b, w, cvx_optval, net, seed, group_flag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global ell_1 ell_2 delta

arc=size( E, 1); % Number of arcs in the graph
x = cell(nodes,1);
temp = cell(nodes, 1);
tempp = cell(nodes, 1);
L_temp = zeros(nodes, 1);
for i=1:nodes
    x{i} = zeros(n, 1);
end
for i=1:nodes
    L_temp(i) = norm(A{i})^2;
end
G = eye(nodes) - G/(max(diag(G))+1);
Gbar = 1/2*(G + eye(nodes));

%GADM SOLUTION
tic

tol1 = 1e-3;
tol2 = 1e-4;    
opt_flag = 1;
prox_counter = 0;
outer_iter = 0;
Relfun = [];
Infeas = [];
Walltime = [];

xhalf = x;
xp = x;
alpha = 2*min(eig(Gbar))/max(L_temp);

if group_flag
    while(opt_flag)
        xpp = xp;
        xp = x;
        % Solve subproblems
        temp_x = G_x(G, x, 1);
        temp_xp = G_x(Gbar, xpp, 1);
        if outer_iter ==0;
            for i = 1:nodes
                prox_counter = prox_counter + 1;
                temp{i} = A{i}*x{i}-b{i};
                xhalf{i} = temp_x{i} - A{i}'*(sign(temp{i}).*min(delta, abs(temp{i})))*alpha;
                for k = 1:nGroups
                    x_temp = sign( xhalf{i}(groups{k}) ).*pos( abs(xhalf{i}(groups{k})) - ell_1/nodes*alpha );
                    x{i}(groups{k}) = x_temp*pos( 1 - ell_2*w(k)*alpha/nodes/norm(x_temp) );
                end
            end
        else
            for i = 1:nodes
                prox_counter = prox_counter + 1;
                temp{i} = A{i}*x{i}-b{i};
                tempp{i} = A{i}*xpp{i} - b{i};
                xhalf{i} = temp_x{i} - temp_xp{i} + xhalf{i} - ( A{i}'*(sign(temp{i}).*min(delta, abs(temp{i})) -...
                    sign(tempp{i}).*min(delta, abs(tempp{i})) ) )*alpha;
                for k = 1:nGroups
                    x_temp = sign( xhalf{i}(groups{k}) ).*pos( abs(xhalf{i}(groups{k})) - ell_1/nodes*alpha );
                    x{i}(groups{k}) = x_temp*pos( 1 - ell_2*w(k)*alpha/nodes/norm(x_temp) );
                end
            end
        end
        alpha = 1*alpha;
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
        xpp = xp;
        xp = x;
        % Solve subproblems
        temp_x = G_x(G, x, 1);
        temp_xp = G_x(Gbar, xpp, 1);
        if outer_iter ==0;
            for i = 1:nodes
                prox_counter = prox_counter + 1;
                temp{i} = A{i}*x{i}-b{i};
                xhalf{i} = temp_x{i} - A{i}'*(sign(temp{i}).*min(delta, abs(temp{i})))*alpha;
                for k = 1:nGroups
                    x_temp = sign(xhalf{i}(groups{i, k})).*pos(abs(xhalf{i}(groups{i, k})) - ell_1/nodes*alpha);
                    x{i}(groups{i, k}) = x_temp*pos(1 - ell_2*w(k)*alpha/nodes/norm(x_temp));
                end
            end
        else
            for i = 1:nodes
                prox_counter = prox_counter + 1;
                temp{i} = A{i}*x{i}-b{i};
                tempp{i} = A{i}*xpp{i} - b{i};
                xhalf{i} = temp_x{i} - temp_xp{i} + xhalf{i} - ( A{i}'*(sign(temp{i}).*min(delta, abs(temp{i})) -...
                    sign(tempp{i}).*min(delta, abs(tempp{i})) ) )*alpha;
                for k = 1:nGroups
                    x_temp = sign( xhalf{i}(groups{i, k}) ).*pos( abs(xhalf{i}(groups{i, k})) - ell_1/nodes*alpha );
                    x{i}(groups{i, k}) = x_temp*pos( 1 - ell_2*w(k)*alpha/nodes/norm(x_temp) );
                end
            end
        end
        alpha = 1*alpha;
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
        row_name = ['PG_EXTRA_tree','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    else
        row_name = ['PG_EXTRA_clique','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    end
else
    if net == 1
        row_name = ['PG_EXTRA_tree','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    else
        row_name = ['PG_EXTRA_clique','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    end
end
output_pgextra = dataset({[relfun, dist, t, outer_iter, prox_counter/nodes], 'Rel', 'ConsensusViolation',...
    'WallTime', 'Iterations', 'ProxOperations'}, 'ObsNames', {row_name});
if isempty(dir('Result.csv'))
    output = output_pgextra;
    export(output, 'File', 'Result.csv', 'Delimiter', ',')
else
    output_old = dataset('File', 'Result.csv', 'Delimiter', ',', 'ReadObsNames', true);
    if ismember(row_name, output_old.Properties.ObsNames)
        output_old(row_name, :) = output_pgextra;
        output_new = output_old;
    else
        output_new = [output_old; output_pgextra];
    end
    export(output_new, 'File', 'Result.csv', 'Delimiter', ',')
end

end

