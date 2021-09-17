function minimum=plotFourCheb(u,a,b,nb)
%Plots the function represented by the Fourier-Chebyshev coefficients in f, 
%[a,b]x[0,2*pi] discretized with step size nb(1) and nb(2).
%Normalization: 
%u(r,psi)=\sum_{n} u_n(r) e^{i*n*psi}, where
%u_n(r)=u_{0,n} + 2\sum_{k\geq 1} u_{k,n} T_k(2*x/(b-a)+(b+a)/(b-a)).

if nargin<4
    nb=[100,100];
end
if nargin<2
    a=-1;
    b=1;
end

[K,N]=size(u);
N=(N+1)/2;

u(2:end,:)=2*u(2:end,:);
u=reshape(u,[K*(2*N-1),1]);

I_N=speye(2*N-1);

r=(a:(b-a)/nb(1):b)';
MK=cos(acos(2*r/(b-a)-(b+a)/(b-a))*(0:K-1));%Very bad way of coding this
MK=kron(I_N,MK);

psi=(0:2*pi/nb(2):2*pi)';
MN=exp(1i*psi*(-N+1:N-1));
MN=kron(MN,speye(length(r)));

[R,Psi]=meshgrid(r,psi);
R=R';
Psi=Psi';

uK=MK*u;
uKN=MN*uK;

uKN_r=real(uKN);
uKN_i=imag(uKN);

uKN_r=reshape(uKN_r,[length(r),length(psi)]);
surf(R,Psi,uKN_r)
xlabel('r')
ylabel('\psi')
minimum=min(min(uKN_r));

if max(abs(uKN_i))>10^-5
    uKN_i=reshape(uKN_i,[length(r),length(psi)]);
    figure
    surf(R,Psi,uKN_i)
    xlabel('r')
    ylabel('\psi')
    title('Imaginary part')
end
    