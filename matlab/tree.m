function [ E, G ] = tree( nodes )

E = zeros(nodes-1, nodes); 
G = zeros(nodes, nodes); 
for i=1:nodes-1
    E(i, 1) = 1;
    E(i, i+1 ) = -1;
    if i==1
        G(1, 1) = nodes-1;
        G(1, 2:end) = -1;
    else
        G(i, i) = 1;
        G(i ,1) = -1;
    end
end
G(nodes, nodes) = 1;
G(nodes, 1 ) = -1;

end

