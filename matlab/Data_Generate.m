function [ n, E, G, A, b, groups, w ] = Data_Generate( nodes, gsize, nGroups, net, seed, ...
    group_flag, edges)

    rng(seed, 'twister')
    overlap=0;  % number of variables overlapped between two consecutive groups
    n=nGroups*(gsize-overlap)+overlap;   % totoal number of variables
    m=floor(n/2/nodes);      % sample size (first dimension of martix A_i) 

    if net==1
        [E, G] = tree(nodes);
    elseif net==2
        [E, G] = clique(nodes);
    elseif net==3
        [E, G] = small_world(nodes, edges, seed);
    end
    
    E = sparse(E);
    G = sparse(G);
    
    A = cell(nodes,1);
    b = cell(nodes,1);
    groups = cell(nGroups,1);
    coe=zeros(n, 1);

    [ Abar, bbar ]=gentoy_group(n/2, nGroups, gsize, overlap);
    w = ones(nGroups, 1);
    for i=1:nodes
        if rand<0.5
            A{i} = Abar( (i-1)*m+1:i*m, :)/(2*sqrt(gsize/100));
            b{i} = bbar( (i-1)*m+1:i*m )/(2*sqrt(gsize/100));
        else
            A{i} = 2*Abar( (i-1)*m+1:i*m, :)/(2*sqrt(gsize/100));
            b{i} = 2*bbar( (i-1)*m+1:i*m )/(2*sqrt(gsize/100));
        end
    end

    if group_flag
        for k=1:nGroups
            groups{k}=[(k-1)*gsize+1:k*gsize]';
            coe(groups{k})=w(k);
        end
    else
        for i=1:nodes
            rand_ind = randperm(n)';
            for k=1:nGroups
                groups{i, k} = rand_ind((k-1)*gsize+1:k*gsize);
                coe(groups{i, k}) = w(i);
            end
        end
    end
end

