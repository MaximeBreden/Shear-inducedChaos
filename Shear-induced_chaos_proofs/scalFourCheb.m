function s=scalFourCheb(U,V,Mat)

%Computation of the L^2 scalar product from the Fourier-Cheybyshev
%representation

if nargin<3
    [K,N]=size(U);
    N=(N+1)/2;
    if N==1
        error('If the inputs are already reshaped as a vector, you need to provide the matrix for the scalar product')
    end
    if exist('intval','file') && isintval(U(1))
        Mat=iscalFourCheb_mat(K,N);
    else
        Mat=scalFourCheb_mat(K,N);
    end
end

if size(U,2)==1
    s=V'*Mat*U;    
else
    s=reshape(V,[my_numel(V),1])'*Mat*reshape(U,[my_numel(U),1]);
end
    