function v=derCheb(u)

% The input u is supposed to be a column. The convention is that
% u = u_0 + 2*\sum u_k T_k. The output v contains the Chebyshev
% coefficients of u'

K=length(u);

v=zeros(K,1);
even=(0:2:K-1)';
odd=(1:2:K-1)';

if exist('intval','file') && isintval(u(1))
    v=derCheb_mat(K)*u;%would also be fine with float, but is a bit slower 
                       %than what is below. Intlab doesn't seem to like 
                       %cumsum, which is why we don't use the formula below
                       %if the input is of intval type
else
    if mod(K,2)
        v(even(1:end-1)+1)=2*cumsum(odd.*u(odd+1),'reverse');
        v(odd+1)=2*cumsum(even(2:end).*u(even(2:end)+1),'reverse');
    else
        v(even+1)=2*cumsum(odd.*u(odd+1),'reverse');
        v(odd(1:end-1)+1)=2*cumsum(even(2:end).*u(even(2:end)+1),'reverse');
    end
end
end