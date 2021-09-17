function [lb_1,lb_2,eta_L1,eta_L2]=rough_lb_nbeig_tD(SDE,Domain,target,K,N)
%%%% Rough lower bound on the number of eigenvalues of tD we need, %%%%
%%%% with close to optimal paramter eta_L1 and eta_L2 %%%%

%Computational parameters
if nargin<4
    K=30;%Chebyshev in r
    N=10;%Fourier in psi
end
if nargin<3
    target=0;
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
I_tens=kron(I_N,I_K);

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
norm_f2=(abs(f(1))+2*sum(abs(f(2:end))))^2;


%g
e0N=zeros(2*N-1,1);
e0N(N)=1;
cos1N=zeros(2*N-1,1);
cos1N(N-1)=1/2;
cos1N(N+1)=1/2;
Mg=2*kron(convomat_exp(b*e0N+sqrt(a^2+b^2)*cos1N),Mr2);
norm_r2g2=(rmax^3*(b+sqrt(a^2+b^2)))^2;


% df=derCheb(f);
% df=kron(transpose(e0N),df);
% dg=2*kron(transpose(-sqrt(a^2+b^2)*sin1N),r2);
% h=df/resc+dg;
% figure
% h0_inf=plotFourCheb(h,rmin,rmax);
% title(['inf h = ',num2str(h0_inf)])

%h0=inf(df/dr+dg/dpsi)=inf h
C1=3*a+2*sqrt(a^2+b^2);
C2=sigma^2/2;
h0=alpha-max(C1*rmax^2+C2/rmax^2,C1*rmin^2+C2/rmin^2);

%h1 = sup h, ||h||_\infty = max(|h0|,|h1|)
C1=3*a-2*sqrt(a^2+b^2);
C2=sigma^2/2;
if C1<=0 
    h1=alpha-C1*rmax^2-C2/rmax^2;
else
    r_extr=(C2/C1)^(1/4);
    if (r_extr >= rmin) && (r_extr <= rmax)
        h1=alpha-C1*r_extr^2-C2/r_extr^2;
    else
        h1=alpha-min(C1*rmax^2+C2/rmax^2,C1*rmin^2+C2/rmin^2);
    end
end
norm_h2=max(abs(h0),abs(h1))^2;


%e
e=kron(transpose(e0N),alpha*e0K-2*a*r2)+sqrt(a^2+b^2)*kron(transpose(sin1N),r2);
Me=kron(I_N,alpha*M0K-2*a*Mr2)+sqrt(a^2+b^2)*kron(convomat_exp(sin1N),Mr2);

figure
plotFourCheb(e);
title('The expansion/contraction rate e')

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
PEVal_L=real(PEVal_L);

[PEVect_Lstar,PEVal_Lstar]=eigs(Trunc_tens*Lstar*Mat_DBC_tens,Trunc_tens*Mat_DBC_tens,1,'SM');
PEVect_Lstar=Mat_DBC_tens*PEVect_Lstar;
PEVect_Lstar=reshape(PEVect_Lstar,[K,2*N-1]);
PEVect_Lstar=symmetrize(PEVect_Lstar);
PEVect_Lstar=PEVect_Lstar/real(sqrt(scalFourCheb(PEVect_Lstar,PEVect_Lstar)));
PEVal_Lstar=real(PEVal_Lstar);

figure
val=plotFourCheb(PEVect_L);
if val<-10*eps
    PEVect_L=-PEVect_L;
    clear gcf
    plotFourCheb(PEVect_L);
end
title(['Principal eigenvector of L, associated eigenvalue = ',num2str(PEVal_L)])
figure
val=plotFourCheb(PEVect_Lstar);
if val<-10*eps
    PEVect_Lstar=-PEVect_Lstar;
    clear gcf
    plotFourCheb(PEVect_Lstar);
end
title(['Principal eigenvector of L^*, associated eigenvalue = ',num2str(PEVal_Lstar)])

Lyapunov_approx=scalFourCheb(reshape(Me*reshape(PEVect_L,[numel(PEVect_L),1]),[K,2*N-1]),PEVect_Lstar);
fprintf('\n\nApproximate value of the first Lyapunov exponent (not rigorous) :\n %f\n\n',Lyapunov_approx)

%% Estimate on the number of eigenvalues of tD needed for the validation of eta

tlambda=PEVal_L;
tU=PEVect_L;
tu=reshape(tU,[numel(tU),1]);
norm_tu2=real(scalFourCheb(tU,tU));

Astar=Lstar-tlambda*I_tens;
Astar_tu=Astar*tu;
Astar_tU=reshape(Astar_tu,[K,2*N-1]);
norm_Astar_tu2=real(scalFourCheb(Astar_tU,Astar_tU));

CV2=norm_f2+norm_r2g2;
eta_S=0.9*norm_tu2;

%Finding an "optimal" eta_L
func=@(eta_L) ( (1-eta_L)/eta_L*CV2-tlambda*sigma^2 ...
    + sqrt( ((1-eta_L)/eta_L*CV2-tlambda*sigma^2)^2-4*(1-eta_L)*sigma^4/4*(tlambda^2+tlambda*h0-1/eta_S*norm_Astar_tu2) ) )...
    /(2*(1-eta_L)*sigma^4/4);

eta_L1=fminbnd(func,0,1)

s2=(1-eta_L1)*sigma^4/4;
s1=(1-eta_L1)/eta_L1*CV2-tlambda*sigma^2;
s0=tlambda^2+tlambda*h0-1/eta_S*norm_Astar_tu2-target;
min_eig_tDelta=(s1+sqrt(s1^2-4*s2*s0))/(2*s2);
lb_1=min_eig_tDelta

%% Estimate on the number of eigenvalues of tD needed for the validation of phi

tlambda=PEVal_Lstar;
tU=PEVect_Lstar;
tu=reshape(tU,[numel(tU),1]);
norm_tu2=real(scalFourCheb(tU,tU));

A=L-tlambda*I_tens;
A_tu=A*tu;
A_tU=reshape(A_tu,[K,2*N-1]);
norm_A_tu2=real(scalFourCheb(A_tU,A_tU));

CV2=norm_f2+norm_r2g2;
eta_S=0.9*norm_tu2;

%Finding an "optimal" eta_L
func=@(eta_L2) ( CV2/eta_L2(1)-tlambda*sigma^2 ...
    + sqrt( (CV2/eta_L2(1)-tlambda*sigma^2)^2-4*(1-eta_L2(1)-eta_L2(2))*sigma^4/4*(tlambda^2+tlambda*h0-norm_A_tu2/eta_S-norm_h2/eta_L2(2)-target) ) )...
    /(2*(1-eta_L2(1)-eta_L2(2))*sigma^4/4);

Aineq=[1,1];
bineq=0.99;

% eta_L2=patternsearch(func,[0.1;0.1],Aineq,bineq,[],[],[0,0],[1,1])
% %%% The following alternative to patternsearch does not require the global
% %%% optimization toolbox
eta_L2=fmincon(func,[0.1;0.1],Aineq,bineq,[],[],[0,0],[1,1]) 

s2=(1-eta_L2(1)-eta_L2(2))*sigma^4/4;
s1=CV2/eta_L2(1)-tlambda*sigma^2;
s0=tlambda^2+tlambda*h0-norm_A_tu2/eta_S-norm_h2/eta_L2(2)-target;
min_eig_tDelta=(s1+sqrt(s1^2-4*s2*s0))/(2*s2);
lb_2=min_eig_tDelta
