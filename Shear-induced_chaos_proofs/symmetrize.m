function v=symmetrize(U)

%Making sure that the function represented by the Fourier-Chebyshev
%coefficients in U is real-valued.

N=size(U,2);
N=(N+1)/2;

v=(U(:,N:end)+conj(flip(U(:,1:N),2)))/2;
v(:,1)=real(v(:,1));
v=[conj(flip(v(:,2:end),2)) v];