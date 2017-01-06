function [X, Y, T, Tw, w] = gentoy_group(n, ng, g_size, overlap)
% X: design matrix
% Y: output
% T: ng times p, indicate each group contains which variables 
% Tw: weight for each group
% w: true regression coefficients
    
    d = ng*(g_size - overlap) + overlap;

    r_T = kron(1:ng, ones(1, g_size));
    c_T = zeros(ng*g_size, 1);
    s = 1;
    for g = 1:ng
        c_T((g-1)*g_size+1:g*g_size, 1) = s:s+g_size-1;
        s = s + g_size - overlap;
    end
    T = sparse(r_T, c_T, 1, ng, d);

    Tw = ones(ng,1);  % uniform weight

    tmp = 1:d;
    w = ((-1).^tmp(:)).*(exp(-(tmp(:)-1)/g_size));  
    sn2 = 0;  % signal to noise ratio

    X = randn(n,d);
    Y = X*w + sn2*rand(n,1);     
end    
