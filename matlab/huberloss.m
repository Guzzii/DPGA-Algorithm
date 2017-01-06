function [ y ] = huberloss( x, delta )

y = zeros(length(x), 1);

ind_1 = find( abs(x)<=delta);
ind_2 = find( abs(x)>delta);
y(ind_1) = 1/2*x(ind_1).^2;
y(ind_2) = delta*( abs(x(ind_2)) - delta/2 );

end

