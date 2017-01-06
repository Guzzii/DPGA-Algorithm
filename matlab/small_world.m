function [ E, G ] = small_world( nodes, edges, seed )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    rng(seed, 'twister');
    
    G = zeros(nodes, nodes);
    E = zeros(nodes+edges, nodes);
    if nodes<=3
        error('# of nodes needs to be >3')
    end
    
    ind = [];
    ind = [ind, nodes-3];
    
    for i =1:nodes-1
        G(i, i+1) = -1;
        if i<=nodes-3;
            ind = [ind, ind(end)+nodes-2-i];
        end
    end
    G(1, nodes) = -1;
     
    loc = randperm((nodes^2-3*nodes)/2);
    loc = loc(1: edges);
    
    for edge=1:edges
        row = find( ind>=loc(edge), 1, 'first');
        if row==1
            col = row+2+loc(edge)-1;
        else
            col = row+2+loc(edge)-ind(row-1)-1;
        end
        
        G(row, col) = -1;
    end
    G = G + G';
    G(1:nodes+1:end) = -sum(G, 2);
    
    k = 0;
    for i = 1:nodes-1
        ind_to = find(G(i, :));
        for j = 1:length(ind_to)
            if ind_to(j)>i
                k = k+1;
                E(k, i) = 1;
                E(k, ind_to(j)) = -1;
            end
        end
    end
end

