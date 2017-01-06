function [ E, G ] = clique( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = zeros(2*N, N); 
G = -ones(N, N); 
k=0;
for i=1:N
    for j=1:N
        if j>i
            k = k+1;
            E(k, i) = 1;
            E(k, j ) = -1;
        end
    end
    G(i, i) = N - 1;
end
end
