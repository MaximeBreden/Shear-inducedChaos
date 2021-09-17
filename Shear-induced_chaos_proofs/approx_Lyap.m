function [Lyap,eta,phi]=approx_Lyap(SDE,Domain,K,N,disp_e)
% Non-rigorous computation of the conditioned Lyapunov exponent

if nargin<5
    disp_e=0;
end

%Computational parameters
if nargin<3
    K=30;%Chebyshev in r
    N=10;%Fourier in psi
end

a=SDE.a;
b=SDE.b;
alpha=SDE.alpha;
sigma=SDE.sigma;
rmin=Domain.rmin;
rmax=Domain.rmax;
resc=Domain.resc;

%% Initialization of all the matrices
I_K=speye(K);
I_N=speye(2*N-1);

%derivatives in psi
MDN=1i*spdiags((-N+1:N-1)',0,2*N-1,2*N-1);
MDN_tens=kron(MDN,I_K);
MD2N=-spdiags((-N+1:N-1)'.^2,0,2*N-1,2*N-1);

%derivatives in r
MDK=derCheb_mat(K);
MDK_tens=kron(I_N,MDK);
MD2K=MDK^2;

%vectors 1 and r in Chebyshev
e0K=zeros(K,1);
e0K(1)=1;
M0K=sparse(convomat_coscos(e0K));
e1K=zeros(K,1);
e1K(2)=1/2;

%r (rescaled) and powers of r
r=(rmax-rmin)/2*e1K+(rmax+rmin)/2*e0K;
Mr=convomat_coscos(r);
r2=Mr*r;
Mr2=Mr*Mr;
r3=Mr*r2;
inv_r=Mr\e0K;
Minv_r=convomat_coscos(inv_r);
inv_r2=Minv_r*inv_r;
Minv_r2=convomat_coscos(inv_r2);
% figure(1)
% plotCheb(inv_r2,rmin,rmax);

% sin(psi)
sin1N=zeros(2*N-1,1);
sin1N(N-1)=1i/2;
sin1N(N+1)=-1i/2;

%the operator \tilde{\Delta} 
tDelta=kron(I_N,MD2K)/resc^2+4*kron(MD2N,Minv_r2);


%f
f=alpha*r-a*r3+sigma^2/2*inv_r;
Mf=kron(I_N,convomat_coscos(f));

%g
e0N=zeros(2*N-1,1);
e0N(N)=1;
cos1N=zeros(2*N-1,1);
cos1N(N-1)=1/2;
cos1N(N+1)=1/2;
Mg=2*kron(convomat_exp(b*e0N+sqrt(a^2+b^2)*cos1N),Mr2);

%e
Me=kron(I_N,alpha*M0K-2*a*Mr2)+sqrt(a^2+b^2)*kron(convomat_exp(sin1N),Mr2);

if disp_e
    e=kron(transpose(e0N),alpha*e0K-2*a*r2)+sqrt(a^2+b^2)*kron(transpose(sin1N),r2);
    figure
    plotFourCheb(e,rmin,rmax);
    title('The expansion/contraction rate e')
end

%V, Vstar
V=Mf/resc*MDK_tens+Mg*MDN_tens;
Vstar=-(MDK_tens*Mf/resc+MDN_tens*Mg);

%L, Lstar, cL
L=sigma^2/2*tDelta+V;
Lstar=sigma^2/2*tDelta+Vstar;

%Dirichlet boundary condition in r
Mat_DBC=[zeros(2,K-2);speye(K-2)];
Mat_DBC(1,1:2:end)=-2;
Mat_DBC(2,2:2:end)=-1;
Mat_DBC_tens=kron(I_N,Mat_DBC);

Trunc=[speye(K-2),zeros(K-2,2)];
Trunc_tens=kron(I_N,Trunc);


%% Getting the principal eigenvalue/eigenvector of L and L*
[PEVect_L,PEVal_L]=eigs(Trunc_tens*L*Mat_DBC_tens,Trunc_tens*Mat_DBC_tens,1,'SM');
PEVect_L=Mat_DBC_tens*PEVect_L;
PEVect_L=reshape(PEVect_L,[K,2*N-1]);
PEVect_L=symmetrize(PEVect_L);
PEVect_L=PEVect_L/real(sqrt(scalFourCheb(PEVect_L,PEVect_L)));
PEVal_L=real(PEVal_L)

[PEVect_Lstar,PEVal_Lstar]=eigs(Trunc_tens*Lstar*Mat_DBC_tens,Trunc_tens*Mat_DBC_tens,1,'SM');
PEVect_Lstar=Mat_DBC_tens*PEVect_Lstar;
PEVect_Lstar=reshape(PEVect_Lstar,[K,2*N-1]);
PEVect_Lstar=symmetrize(PEVect_Lstar);
PEVect_Lstar=PEVect_Lstar/real(sqrt(scalFourCheb(PEVect_Lstar,PEVect_Lstar)));
PEVal_Lstar=real(PEVal_Lstar)

figure
val=plotFourCheb(PEVect_L);
if val<-0.1
    PEVect_L=-PEVect_L;
end
close gcf
figure
val=plotFourCheb(PEVect_Lstar);
if val<-0.1
    PEVect_Lstar=-PEVect_Lstar;
end
close gcf

% if PEVal_L>=0 || PEVal_Lstar>=0 || abs((PEVal_L-PEVal_Lstar)/PEVal_L)>10^-4
%     figure
%     plotFourCheb(PEVect_L,rmin,rmax);
%     title(['\lambda_{\eta} = ',num2str(PEVal_L)])
%     figure
%     plotFourCheb(PEVect_Lstar,rmin,rmax);
%     title(['\lambda_{\phi} = ',num2str(PEVal_Lstar)])
%     error('Something might have gone wrong with the numerical computation. Maybe try increasing K and/or N')
% end

Lyap=real(scalFourCheb(reshape(Me*reshape(PEVect_L,[numel(PEVect_L),1]),[K,2*N-1]),PEVect_Lstar))
eta=PEVect_L;
phi=PEVect_Lstar;