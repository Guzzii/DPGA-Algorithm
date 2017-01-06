function [ out, cvx_optval ] = prox1( xp, p, y, G, j, opt )
% Proximal operation for regularizer and huber function.
%   Detailed explanation goes here
global n A b groups w ell_1 ell_2 c nodes nGroups 

yj = zeros(n, 1);

if opt==1
    obj='';
    for k=1:nGroups
        obj = [obj,'ell_2/nodes*w(',num2str(k),')*norm(xc(groups{',num2str(k),'}))+'];
    end
    obj = [obj, 'ell_1/nodes*norm(xc,1)+0.5*sum(huber(A{',num2str(j),'}*xc-b{',num2str(j),'}))' ];

    indList = find( G(:, j)~=0 );
    aj = sum( G(indList, j).^2 );
    for i=1:length(indList)
        yj =yj - ( G(indList(i), j)*y{indList(i)} - G(indList(i), j)^2*xp ) - G(indList(i), j)*p{indList(i)}/c;
    end
    obj = [obj, '+0.5*c*aj*sum_square( xc-yj/aj )'];

    cvx_begin quiet
        variable xc(n)
        minimize(obj)
        subject to
    cvx_end

    out = xc;
    cvx_optval;
    
elseif opt==2
    obj='';
    for k=1:nGroups
        obj = [obj,'ell_2/nodes*w(',num2str(k),')*norm(xc(groups{',num2str(j),',',num2str(k),'}))+'];
    end
    obj = [obj, 'ell_1/nodes*norm(xc,1)+0.5*sum(huber(A{',num2str(j),'}*xc-b{',num2str(j),'}))' ];

    indList = find( G(:, j)~=0 );
    aj = sum( G(indList, j).^2 );
    for i=1:length(indList)
        yj =yj - ( G(indList(i), j)*y{indList(i)} - G(indList(i), j)^2*xp ) - G(indList(i), j)*p{indList(i)}/c;
    end
    obj = [obj, '+0.5*c*aj*sum_square( xc-yj/aj )'];

    cvx_begin quiet
        variable xc(n)
        minimize(obj)
        subject to
    cvx_end

    out = xc;
    cvx_optval;
    
end

