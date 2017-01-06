function [ s ] = G_x( G, x, opt )
%   Computes s iterates of DPGA-I and DPGA-II
global n

k = size(G, 1);
indList = cell(k, 1);
s = cell(k, 1);

for i=1:k
    indList{i} = find(G(i, :)~=0);
    s{i} = zeros(n, 1);
end

if opt==1
    for i=1:k
        for j=1:length(indList{i})
            s{i} = s{i} + G(i, indList{i}(j))* x{indList{i}(j)} ;
        end
    end
    
elseif opt==2
    for i=1:k
        for j=1:length(indList{i})
            s{i} = s{i} + G(i, indList{i}(j))* x{indList{i}(j) } /(G(i, i)+1) ;
        end
    end
end

end

