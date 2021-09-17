function operator_S2=initialize_S2(SDE,Domain,tU,tlambda,op_tD,op_L,op_S1,approx_r,weights,K,N)
% Matrix representations of the operators remaining operators involved in SS^* (for L^*)

a=SDE.a;
b=SDE.b;
alpha=SDE.alpha;
sigma=SDE.sigma;
ia=SDE.ia;
ib=SDE.ib;
ialpha=SDE.ialpha;
isigma=SDE.isigma;

rmin=Domain.rmin;
rmax=Domain.rmax;
irmin=Domain.irmin;
irmax=Domain.irmax;
resc=Domain.resc;
iresc=Domain.iresc;

I_K=speye(K);
I_N=speye(2*N-1);
I_tens=kron(I_N,I_K);
I_tens_ext=speye((3*K-2)*(2*N+3));

tu=reshape(tU,[numel(tU),1]);
norm_tu2=real(scalFourCheb(tU,tU));

itlambda=intval(tlambda);
itu=intval(tu);
itu=op_L.iMat_BC_tens*op_L.Reset_for_BC_tens*itu;
itU=reshape(itu,[K,2*N-1]);
itU=symmetrize(itU);
itU=[[intval(zeros(K,2)),itU,intval(zeros(K,2))];intval(zeros(2*(K-1),2*N+3))];
inorm_tu2=real(scalFourCheb(itU,itU));
itu=reshape(itU,[(3*K-2)*(2*N+3),1]);

A=op_L.L-tlambda*I_tens;
Astar=op_L.Lstar-tlambda*I_tens;
A_tu=A*tu;
A_tU=reshape(A_tu,[K,2*N-1]);
norm_A_tu2=real(scalFourCheb(A_tU,A_tU));

iA=op_L.iL-itlambda*I_tens_ext;
iAstar=op_L.iLstar-itlambda*I_tens_ext;
iA_tu=iA*itu;
iA_tU=reshape(iA_tu,[3*K-2,2*N+3]);
inorm_A_tu2=real(scalFourCheb(iA_tU,iA_tU));

xi0=weights.xi0;
ixi0=weights.ixi0;
eta_L2=weights.eta_L2;
ieta_L2=intval(eta_L2);

CV=sqrt(op_L.norm_f2+op_L.norm_r2g2);
iCV=sqrt(op_L.inorm_f2+op_L.inorm_r2g2);
eta_S=0.9*norm_tu2;
ieta_S=intval(eta_S);

s2=(1-eta_L2(1)-eta_L2(2))*sigma^4/4;
s1=CV^2/eta_L2(1)-tlambda*sigma^2;
s0=tlambda^2+tlambda*op_L.h0-norm_A_tu2/eta_S-op_L.norm_h2/eta_L2(2);
slambda=xi0*(norm_tu2-eta_S);

is2=(1-ieta_L2(1)-ieta_L2(2))*isigma^4/4;
is1=iCV^2/ieta_L2(1)-itlambda*isigma^2;
is0=itlambda^2+itlambda*op_L.ih0-inorm_A_tu2/ieta_S-op_L.inorm_h2/ieta_L2(2);
islambda=ixi0*(inorm_tu2-ieta_S);

%||h + lambda||_\infty 
inorm_hplambda=max(abs(op_L.ih0+itlambda),abs(op_L.ih1+itlambda));

coef_C0_tens=repmat([1 2*ones(1,3*(K-1))],[1,2*N+3]);   

%%Boundary condition in r
%For the end problem    
frmax=alpha*rmax-a*rmax^3+sigma^2/(2*rmax);
frmin=alpha*rmin-a*rmin^3+sigma^2/(2*rmin);
Mat_DBC2=[zeros(4,K-4);speye(K-4)];
DBC2=[ones(1,K);
      (-1).^(0:K-1);
      sigma^2/2*(0:K-1).^2.*((0:K-1).^2-1)/3/resc^2-frmax*(0:K-1).^2/resc;
      (-1).^(0:K-1).*(sigma^2/2.*(0:K-1).^2.*((0:K-1).^2-1)/3/resc^2+frmin*(0:K-1).^2/resc)]; 
DBC2(:,2:end)=2*DBC2(:,2:end);

DBC2_to_inv=DBC2(:,1:4);
DBC2_remainder=DBC2(:,5:end);            
Mat_DBC2(1:4,:)=-DBC2_to_inv\DBC2_remainder;
Mat_DBC2_tens=kron(I_N,Mat_DBC2);
Mat_DBC2_tens=[Mat_DBC2_tens zeros(K*(2*N-1),1);
               zeros(1,(K-4)*(2*N-1)) 1];

ifrmax=ialpha*irmax-ia*irmax^3+isigma^2/(2*irmax);
ifrmin=ialpha*irmin-ia*irmin^3+isigma^2/(2*irmin);
iMat_DBC2=intval([zeros(4,K-4);speye(K-4)]);
iK2=intval(0:K-1).^2;
iDBC2=[ones(1,K);
      (-1).^(0:K-1);
      isigma^2/2*iK2.*(iK2-1)/3/iresc^2-ifrmax*iK2/resc;
      (-1).^(0:K-1).*(isigma^2/2.*iK2.*(iK2-1)/3/iresc^2+ifrmin*iK2/iresc)]; 
iDBC2(:,2:end)=2*iDBC2(:,2:end);
iDBC2_to_inv=iDBC2(:,1:4);
iDBC2_remainder=iDBC2(:,5:end);            
iMat_DBC2(1:4,:)=-iDBC2_to_inv\iDBC2_remainder;
iMat_DBC2_tens=my_kron(I_N,iMat_DBC2);
iMat_DBC2_tens=[iMat_DBC2_tens zeros(K*(2*N-1),1);
                zeros(1,(K-4)*(2*N-1)) 1];     

%% Preparations for the error terms
I_K_ext=speye(3*K-2);
I_N_ext=speye(2*N+3);

%derivatives in psi
MDN_ext=1i*spdiags((-N-1:N+1)',0,2*N+3,2*N+3);

%derivatives in r
MDK_ext=derCheb_mat(3*K-2);

%full der on u
iMDN_ext=intval(MDN_ext);
iMDK_ext=intval(MDK_ext);

iD4Psi=my_kron(iMDN_ext^4,I_K_ext);
iD3Psi=my_kron(iMDN_ext^3,I_K_ext);
iD2Psi=my_kron(iMDN_ext^2,I_K_ext);
iDPsi=my_kron(iMDN_ext,I_K_ext);
iD3r=my_kron(I_N_ext,iMDK_ext^3/iresc^3);
iD2r=my_kron(I_N_ext,iMDK_ext^2/iresc^2);
iDr=my_kron(I_N_ext,iMDK_ext/iresc);
iD2rD2Psi=my_kron(iMDN_ext^2,iMDK_ext^2/iresc^2);
iDrD2Psi=my_kron(iMDN_ext^2,iMDK_ext/iresc);
iDrDPsi=my_kron(iMDN_ext,iMDK_ext/iresc);

%r stuff
% iMinv_r=convomat_coscos(iinv_r,3*K-2);
% iMinv_r_tens=my_kron(I_N_ext,iMinv_r);
iMinv_r2=convomat_coscos(approx_r.iinv_r2,3*K-2);
iMinv_r2_tens=my_kron(I_N_ext,iMinv_r2);

%first and zero order terms
iMf=op_L.iMf;
iMg=op_L.iMg;
iMh=op_L.iMh;

%% Structure
operator_S2=op_S1;
clear op_S1

operator_S2.s2=s2;
operator_S2.s1=s1;
operator_S2.s0=s0;
operator_S2.slambda=slambda;

operator_S2.is2=is2;
operator_S2.is1=is1;
operator_S2.is0=is0;
operator_S2.islambda=islambda;

tDelta=-op_tD.H;
operator_S2.H0=[s2*tDelta^2+s1*tDelta+s0*I_tens,zeros(K*(2*N-1),1);
              zeros(1,K*(2*N-1)),slambda];

itDelta=op_L.itDelta;
operator_S2.iH0=[is2*itDelta^2+is1*itDelta+is0*I_tens_ext,intval(zeros((3*K-2)*(2*N+3),1));
               intval(zeros(1,(3*K-2)*(2*N+3))),islambda];

operator_S2.H=[A*Astar+xi0*(tu*tu')*scalFourCheb_mat(K,N),-A_tu;
             -xi0*A_tu'*scalFourCheb_mat(K,N),xi0*norm_tu2];

operator_S2.iH=[iA*iAstar+ixi0*(itu*itu')*iscalFourCheb_mat(3*K-2,N+2),-iA_tu;
              -ixi0*iA_tu'*iscalFourCheb_mat(3*K-2,N+2),ixi0*inorm_tu2];

operator_S2.Mat_BC_tens=Mat_DBC2_tens;
operator_S2.iMat_BC_tens=iMat_DBC2_tens;

%stuff needed for the NK theorem
operator_S2.iCV=iCV;
operator_S2.itlambda=itlambda;
operator_S2.itU=itU;
operator_S2.itu=itu;
operator_S2.inorm_tu2=inorm_tu2;
operator_S2.iD2Psi=iD2Psi;
operator_S2.iDr=iDr;
operator_S2.inorm_hplambda=inorm_hplambda;
operator_S2.I_tens_ext=I_tens_ext;
           
%Coefficients for the error terms in A0, A1 and A2
operator_S2.icoef_S0_inv_r2=is2*(4*iD2rD2Psi+16*iMinv_r2_tens*iD4Psi+4*iD2Psi*itDelta)...
                   +is1*(4*iD2Psi);
operator_S2.icoef_S0_der_inv_r2=is2*8*iDrD2Psi;
operator_S2.icoef_S0_der2_inv_r2=is2*4*iD2Psi;
operator_S2.icoef_S0_inv_r2_inv_r2=is2*16*iD4Psi;

operator_S2.icoef_cA_inv_r2=isigma^4*iD2rD2Psi+4*isigma^4*iMinv_r2_tens*iD4Psi+2*isigma^2*iMf*iDrD2Psi...
                            +2*isigma^2*iMg*iD3Psi-2*isigma^2*(itlambda*I_tens_ext+isigma^2*iMinv_r2_tens)*iD2Psi...
                            -isigma^4/4*iD2r-isigma^2/2*iMf*iDr-isigma^2/2*iMg*iDPsi-isigma^2/2*I_tens_ext...
                            +2*isigma^2*iD2Psi*iAstar;
operator_S2.icoef_cA_der_inv_r2=2*isigma^4*iDrD2Psi-isigma^4/4*iDr+2*isigma^2*iMf*iD2Psi-isigma^2/2*iMf;
operator_S2.icoef_cA_der2_inv_r2=isigma^4*iD2Psi;
operator_S2.icoef_cA_inv_r=-isigma^4/4*iD3r-isigma^4*iMinv_r2_tens*iDrD2Psi-isigma^2/2*iMf*iD2r...
               -isigma^2/2*iMg*iDrDPsi+isigma^2/2*itlambda*iDr...
               +isigma^2/2*iDr*iAstar;
operator_S2.icoef_cA_der_inv_r=-isigma^4/2*iD2r-isigma^2/2*iMf*iDr;
operator_S2.icoef_cA_der2_inv_r=-isigma^4/4*iDr;

operator_S2.icoef_cA_inv_r2_inv_r2=4*isigma^4*iD4Psi-isigma^4*iD2Psi;
operator_S2.icoef_cA_inv_r2_inv_r=-isigma^4/4*iDr;
operator_S2.icoef_cA_inv_r_inv_r=-isigma^4/4*iD2r;
operator_S2.icoef_cA_inv_r_der_inv_r2=isigma^4*iD2Psi-isigma^4/4*I_tens_ext;
operator_S2.icoef_cA_inv_r_der_inv_r=-isigma^4/4*iDr;

operator_S2.icoef_A_tu_inv_r2_fullC0=coef_C0_tens*abs(2*isigma^2*iD2Psi*itu);
operator_S2.icoef_A_tu_inv_r_fullC0=coef_C0_tens*abs(isigma^2/2*iDr*itu);           
