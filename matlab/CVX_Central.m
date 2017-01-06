function [ cvx_optval, xc ] = CVX_Central( nodes, gsize, n, nGroups, groups, E, G, A, b, w, seed, ...
    group_flag, alg_flag )

global ell_1 ell_2 delta

%CVX Solution
if group_flag %solvig for same grouping
    obj='';
    for k=1:nGroups
        obj=[obj,'ell_2*w(',num2str(k),')*norm(xc(groups{',num2str(k),'}))+'];
    end
    obj=[obj,'ell_1*norm(xc,1)'] ;
    for i=1:nodes
        obj = [obj, '+0.5*sum(huber(A{',num2str(i),'}*xc-b{',num2str(i),'}))' ];
    end
else %solvig for different grouping
    obj='';
    for i=1:nodes
        for k=1:nGroups
            obj=[obj,'ell_2/nodes*w(',num2str(k),')*norm(xc(groups{',num2str(i),',',num2str(k),'}))+'];
        end
    end
    obj=[obj,'ell_1*norm(xc,1)'] ;
    for i=1:nodes
        obj = [obj,'+0.5*sum(huber(A{',num2str(i),'}*xc-b{',num2str(i),'}))' ];
    end
end

tic
cvx_begin quiet
if alg_flag
    cvx_solver sdpt3
else
    cvx_solver mosek
end
    variable xc(n)
    minimize(obj)
cvx_end
t = toc

if alg_flag
    if group_flag
        row_name = ['SDPT3','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    else
        row_name = ['SDPT3','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    end
else
    if group_flag
        row_name = ['Mosek','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed)];
    else
        row_name = ['Mosek','_',num2str(nodes),'_',num2str(gsize),'_',num2str(seed), '_diff'];
    end
end
output_cvx = dataset({[0, 0, t, NaN, NaN], 'Rel', 'ConsensusViolation', 'WallTime',...
    'Iterations', 'ProxOperations'}, 'ObsNames', {row_name});
if isempty(dir('Result.csv'))
    output = output_cvx;
    export(output, 'File', 'Result.csv', 'Delimiter', ',')
else
    output_old = dataset('File', 'Result.csv', 'Delimiter', ',', 'ReadObsNames', true);
    if ismember(row_name, output_old.Properties.ObsNames)
        output_old(row_name, :) = output_cvx;
        output_new = output_old;
    else
        output_new = [output_old; output_cvx];
    end
    export(output_new, 'File', 'Result.csv', 'Delimiter', ',')
end

end

