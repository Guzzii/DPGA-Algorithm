function [ out, cvx_optval ] = prox2( xp, yp, r, p, s, G, j, opt )
% Proximal operation for regularizer and huber function.
%   Detailed explanation goes here
global n A b groups w ell_1 ell_2 c nodes nGroups 

yj = zeros(n, 1);
out = zeros(n, 1);

if opt==1 %prox for x
    obj = '';
    for k=1:nGroups
        obj = [obj,'ell_2/nodes*w(',num2str(k),')*norm(xc(groups{',num2str(k),'}))+'];
    end
    obj = [obj, 'ell_1/nodes*norm(xc,1)' ];
    indList = find( G(:, j)~=0 );
    aj = sum( G(indList, j).^2 )+1;
    for i=1:length(indList)
        yj =yj - ( G(indList(i), j)*s{indList(i)} - G(indList(i), j)^2*xp ) - G(indList(i), j)*p{indList(i)}/c;
    end
    yj = yj - r + (xp+yp)/2;
    obj = [obj, '+0.5*c*aj*sum_square( xc-yj/aj )'];
    
    cvx_begin quiet
        cvx_solver mosek
        variable xc(n)
        minimize(obj)
        subject to
    cvx_end
    out = xc;
    cvx_optval;
    
elseif opt==2 %prox for x
    obj = '';
    obj = [obj, '0.5*sum(huber(A{',num2str(j),'}*yc-b{',num2str(j),'}))'];
    obj = [obj, '+0.5*c*sum_square( yc - ( r+(xp+yp)/2 ) )'];
    
    cvx_begin quiet
        variable yc(n)
        minimize(obj)
        subject to
    cvx_end
    out = yc;
    cvx_optval;
    
elseif opt==3 %closed form prox for x
    indList = find( G(:, j)~=0 );
    aj = sum( G(indList, j).^2 )+1;
    for i=1:length(indList)
        yj =yj - ( G(indList(i), j)*s{indList(i)} - G(indList(i), j)^2*xp ) - G(indList(i), j)*p{indList(i)}/c;
    end
    yj = (yj - r + (xp+yp)/2)/aj;
    
    for k=1:nGroups
        x_temp = sign( yj(groups{k}) ).*pos( abs(yj(groups{k})) - ell_1/nodes/(c*aj) );
        out(groups{k}) = x_temp*pos( 1 - ell_2*w(k)/nodes/(c*aj)/norm(x_temp) );
    end
    
elseif opt==4 %prox for x in ADMM3
    obj = '';
    obj = [obj, '0.5*sum(huber(A{',num2str(j),'}*yc-b{',num2str(j),'}))'];
    indList = find( G(:, j)~=0 );
    aj = sum( G(indList, j).^2 )+1;
    for i=1:length(indList)
        yj =yj - ( G(indList(i), j)*s{indList(i)} - G(indList(i), j)^2*yp ) - G(indList(i), j)*p{indList(i)}/c;
    end
    yj = yj + r + (xp+yp)/2;
    obj = [obj, '+0.5*c*aj*sum_square( yc-yj/aj )'];
    
    cvx_begin quiet
        variable yc(n)
        minimize(obj)
        subject to
    cvx_end
    out = yc;
    cvx_optval;

elseif opt==5 %closed form prox for x with different groups
    indList = find( G(:, j)~=0 );
    aj = sum( G(indList, j).^2 )+1;
    for i=1:length(indList)
        yj =yj - ( G(indList(i), j)*s{indList(i)} - G(indList(i), j)^2*xp ) - G(indList(i), j)*p{indList(i)}/c;
    end
    yj = (yj - r + (xp+yp)/2)/aj;
    
    for k=1:nGroups
        x_temp = sign( yj(groups{j, k}) ).*pos( abs(yj(groups{j, k})) - ell_1/nodes/(c*aj) );
        out(groups{j, k}) = x_temp*pos( 1 - ell_2*w(k)/nodes/(c*aj)/norm(x_temp) );
    end
end

end

