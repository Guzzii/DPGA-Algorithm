function[x, Relfun, Infeas, Walltime, normsqd, L] = DPGA(nodes, gsize, n, ...
    nGroups, groups, E, G, A, b, w, cvx_optval, xc, net, seed, ...
    group_flag, timlim, eta)
% Solving distributed problem using DPGA-II ADMM method

global ell_1 ell_2 delta

if ~exist('eta', 'var')
   eta = 1; 
end

arc=size(E, 1); % Number of arcs in the graph

c = 3.5156*sqrt(4/arc/min(diag(G))) * ones(nodes, 1);

x = cell(nodes,1);
p = cell(nodes, 1);
s = cell(nodes, 1); 
grad=cell(nodes, 1);
temp = cell(nodes, 1);
subgrad = zeros(n, 1);
f_val = zeros(nodes, 1);
L_true = zeros(nodes, 1);
for i=1:nodes
    x{i} = zeros(n, 1);
    p{i} = zeros(n, 1);
    L_true(i) = norm(A{i})^2;
end
G_sp = kron(G, speye(n));

tic

tol1 = 1e-3;
tol2 = 1e-4;    
opt_flag = 1;
outer_iter = 0;
prox_count = 0;

Relfun = [];
Infeas = [];
Walltime = [];

if eta == 1
    disp('Constant stepsize Strategy');
    backtrack = false;
    L_temp = L_true;
elseif eta > 1
    disp('Adaptive stepsize Strategy')
    backtrack = true;
    L_temp = L_true/eta;
else
    error('eta should be >= 1')
end

L = L_temp+c.*diag(G);

while(opt_flag)
    xp = x;
    if outer_iter==0
        s = mat2cell(G_sp*cell2mat(x), repmat(n, nodes, 1));
    end

    for i=1:nodes
        f_valp = sum( huberloss( A{i}*xp{i}-b{i}, delta ) );

        temp{i} = A{i}*xp{i}-b{i};
        grad{i} = A{i}'*(sign(temp{i}).*min(delta, abs(temp{i})));
        temp{i} = [];
        
        % back-tracking Lipchitz constant
        backtrack_flag = 1;
        while(backtrack_flag)
            x{i} = xp{i} - (grad{i} + c(i)/2*(p{i}+s{i}))/L(i);
            prox_count = prox_count + 1;
            for k=1:nGroups
                if group_flag
                    x_temp = sign(x{i}(groups{k})).*max(abs(x{i}(groups{k}))...
                        - ell_1/nodes/L(i), 0);
                    x{i}(groups{k}) = x_temp*max(1 ...
                        - ell_2*w(k)/nodes/L(i)/norm(x_temp), 0);
                else
                    x_temp = sign(x{i}(groups{i, k})).*max(...
                        abs(x{i}(groups{i, k})) - ell_1/nodes/L(i), 0);
                    x{i}(groups{i, k}) = x_temp*max(1 - ell_2*w(k)...
                        /nodes/L(i)/norm(x_temp), 0);
                end                    
            end
            f_val(i) = sum(huberloss(A{i}*x{i}-b{i}, delta));
            if (backtrack == false) || (L_temp(i) >= L_true(i))
                backtrack_flag = 0;
            else
                % evaluate back-tracking conditioin
                diff = x{i}-xp{i};
                grad_diff = grad{i}'*diff;
                norm_sq_diff = diff'*diff;
                if f_val(i) <= f_valp + grad_diff+0.5*L_temp(i)*norm_sq_diff
                    L_temp(i) = L_temp(i)/eta;
                    backtrack_flag = 0;
                else
                    L_temp(i) = min(L_temp(i)*eta, L_true(i));
                end
                L(i) = L_temp(i)+c(i)*G(i,i);
            end
        end
    end

    % dual update
    s = mat2cell(G_sp*cell2mat(x), repmat(n, nodes, 1));
    for i=1:nodes
        p{i} = p{i}+s{i};
    end

    % compute objective function value
    val1 = 0;
    for i=1:nodes
        for k=1:nGroups
            if group_flag
                val1 = val1+ell_2/nodes*w(k)*norm(x{i}(groups{k}));
            else
                val1 = val1+ell_2/nodes*w(k)*norm(x{i}(groups{i, k}));
            end
        end
        val1 = val1+ell_1/nodes*norm(x{i},1) + f_val(i);
    end
    relfun = abs(cvx_optval - val1)/cvx_optval;

    % compute infeasiblity
    dist=0;
    for i=1:arc
        arcList=find(E(i,:)~=0);
        dist=max( dist, norm(x{arcList(1)}-x{arcList(2)}, 2));
    end
    dist = dist/sqrt(n);

    if relfun < tol1 && dist<tol2
        opt_flag = 0;
    end
    
    t_break = toc;
    
    Relfun = [Relfun, relfun];
    Infeas = [Infeas, dist];
    Walltime = [Walltime, t_break];

    if t_break>timlim
        opt_flag = 0;
    end

    outer_iter = outer_iter+1;

    if mod(outer_iter, 50)==0
        disp(['Out_Ite: ', num2str(outer_iter), ' -- Obj: ', ...
            num2str( val1 ), ' -- |xi=xj|: ', num2str(dist)]);
    end
end

t = toc;

% compute subgradient
normsqd = 0;
for i = 1:nodes
    temp{i} = A{i}*xc-b{i};
    grad{i} = A{i}'*(sign(temp{i}).*min(delta, abs(temp{i})));
    for k=1:nGroups
        if group_flag
            if x{i}(groups{k})/sqrt(length(groups{k}))<=1e-6
                tempp = pos(abs(grad{i}(groups{k})) - ell_1/nodes);
                subgrad(groups{k}) = sign(grad{i}(groups{k})).*tempp*max(0,...
                    1-lambda*ell_2*w(k)/(nodes*norm(tempp,2)));
            else
                tempp =  grad{i}(groups{k}) + ell_2*w(k)*x{i}(groups{k})...
                    /(nodes*norm(x{i}(groups{k}),2));
                ind_z = abs(x{i}(groups{k})) <= 1e-6;
                ind_nz = abs(x{i}(groups{k})) > 1e-6;
                subgrad(groups{k}) = (tempp + ell_1/nodes...
                    *sign(x{i}(groups{k}))).*ind_nz ...
                    + sign(tempp).*pos(abs(tempp) - ell_1/nodes ).*ind_z;
            end
        else
            if x{i}(groups{i, k})/sqrt(length(groups{i, k}))<=1e-6
                tempp = pos(abs(grad{i}(groups{i, k})) - ell_1/nodes);
                subgrad(groups{i, k}) = sign(grad{i}(groups{i, k}))....
                    *tempp*max(0,1-lambda*ell_2*w(k)/(nodes*norm(tempp,2)));
            else
                tempp =  grad{i}(groups{i, k}) + ell_2*w(k)*...
                    x{i}(groups{i, k})/(nodes*norm(x{i}(groups{i, k}),2));
                ind_z = abs(x{i}(groups{i, k})) <= 1e-6;
                ind_nz = abs(x{i}(groups{i, k})) > 1e-6;
                subgrad(groups{i, k}) = (tempp ...
                    + ell_1/nodes*sign(x{i}(groups{i, k}) ) ).*ind_nz...
                    +sign(tempp).*pos(abs(tempp) - ell_1/nodes ).*ind_z;
            end
        end
    end
    normsqd = normsqd + subgrad'*subgrad;
end


if net == 1
    if backtrack
        row_name = ['DPGA_II_AS_tree','_',num2str(nodes),'_',...
            num2str(gsize),'_',num2str(seed)];
    else
        row_name = ['DPGA_II_CS_tree','_',num2str(nodes),'_',...
            num2str(gsize),'_',num2str(seed)];
    end
else
    if backtrack
        row_name=['DPGA_II_AS_clique','_',num2str(nodes),'_',...
            num2str(gsize),'_',num2str(seed)];
    else
        row_name=['DPGA_II_CS_clique','_',num2str(nodes),'_',...
            num2str(gsize),'_',num2str(seed)];
    end
end
if group_flag == 0
    row_name = [row_name,'_diff']; 
end

output_dpga_ii = dataset({[relfun, dist, t, outer_iter, prox_count/nodes],...
    'Rel', 'ConsensusViolation', 'WallTime', 'Iterations', 'ProxOperations'},...
    'ObsNames', {row_name});
if isempty(dir('Result.csv'))
    output = output_dpga_ii;
    export(output, 'File', 'Result.csv', 'Delimiter', ',')
else
    output_old = dataset('File', 'Result.csv', 'Delimiter', ',',...
       'ReadObsNames', true);
    if ismember(row_name, output_old.Properties.ObsNames)
        output_old(row_name, :) = output_dpga_ii;
        output_new = output_old;
    else
        output_new = [output_old; output_dpga_ii];
    end
    export(output_new, 'File', 'Result.csv', 'Delimiter', ',')
end

end

