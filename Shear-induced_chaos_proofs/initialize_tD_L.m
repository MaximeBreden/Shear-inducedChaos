function [operator_tD,operator_L,approx_r]=initialize_tD_L(SDE,Domain,K,N)
% All the matrix representations of the operators involved in
% \tilde{\Delta}, L, L*, and the error terms for the Chebyshev series
% enclosure of 1/r and 1/r^2

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

%derivatives in psi
MDN=1i*spdiags((-N+1:N-1)',0,2*N-1,2*N-1);
MDN_tens=kron(MDN,I_K);
MD2N=-spdiags((-N+1:N-1)'.^2,0,2*N-1,2*N-1);

%derivatives in r (on [-1,1], need to add a factor 1/resc for each r
%derivative when on [r_min,r_max] !!!)
MDK=derCheb_mat(K);
MDK_tens=kron(I_N,MDK);
MD2K=MDK^2;

%vectors 1 and r in Chebyshev (on [-1,1])
e0K=zeros(K,1);
e0K(1)=1;
e1K=zeros(K,1);
e1K(2)=1/2;

%rescaled r and powers of r (on [rmin,rmax])
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

%A posteriori validation of inv_r and inv_r2 (see Appendix F of the paper)
inumax=(irmax+irmin)/(irmax-irmin)*(1+sqrt(1-((irmax-irmin)/(irmax+irmin))^2));
inu_r=0.8*1+0.2*inumax %We made no effort to find an "optimal" nu, anything between 1 and numax could be tried.
ir=(irmax-irmin)/2*e1K+(irmax+irmin)/2*e0K;
iMr=convomat_coscos(ir);
iinv_r=intval(inv_r);
contraction_factor=[e0K;0;0]-convo_coscos(ir,iinv_r,K+2);
contraction_norm=[1 2*inu_r.^(1:K+1)]*abs(contraction_factor);
if sup(contraction_norm)<1
    residual=convo_coscos(contraction_factor,iinv_r,2*K+2);
    residual_norm=[1 2*inu_r.^(1:2*K+1)]*abs(residual);
    inv_r_err=sup(residual_norm/(1-contraction_norm));
    disp(['The approximation of 1/r on [r_min,r_max] has been validated with an error bound of ',num2str(inv_r_err)])    
else
    error(['Unable to validate the approximation of 1/r on [r_min,r_max], the contraction factor (which should be less than 1) is equal to ',num2str(contraction_norm.sup),'. Try to decrease inu_r and/or to increase K'])
end

ir2=iMr*ir;
iinv_r2=intval(inv_r2);
contraction_factor2=[e0K;0;0]-convo_coscos(ir2,iinv_r2,K+2);
contraction_norm2=[1 2*inu_r.^(1:K+1)]*abs(contraction_factor2);
if sup(contraction_norm2)<1
    residual2=convo_coscos(contraction_factor2,iinv_r2,2*K+2);
    residual_norm2=[1 2*inu_r.^(1:2*K+1)]*abs(residual2);
    inv_r2_err=sup(residual_norm2/(1-contraction_norm2));
    disp(['The approximation of 1/r^2 on [r_min,r_max] has been validated with an error bound of ',num2str(inv_r2_err)])    
else
    error(['Unable to validate the approximation of 1/r^2 on [r_min,r_max], the contraction factor (which should be less than 1) is equal to ',num2str(contraction_norm2.sup),'. Try to decrease inu_r and/or to increase K'])
end

cnst_der1=cnst_der(inu_r,1)%|| u' ||_{C^0([-1,1])} <= cnst_der1 * || u ||_{\ell^1_{\nu_r}}
cnst_der2=cnst_der(inu_r,(1+inu_r)/2)*cnst_der((1+inu_r)/2,1)%|| u'' ||_{C^0([-1,1])} <= cnst_der2 * || u ||_{\ell^1_{\nu_r}}

%Error estimates
approx_r.iinv_r2=iinv_r2;
approx_r.inv_r_err=inv_r_err;
approx_r.inv_r2_err=inv_r2_err;
approx_r.der_inv_r2_err=cnst_der1/iresc*inv_r2_err;
approx_r.der2_inv_r2_err=cnst_der2/iresc^2*inv_r2_err;
approx_r.der_inv_r_err=cnst_der1/iresc*inv_r_err;
approx_r.der2_inv_r_err=cnst_der2/iresc^2*inv_r_err;

%the operators \tilde{\Delta} and the associated base problem
%\tilde{\Delta}^{(0)}. 
tDelta=kron(I_N,MD2K)/resc^2+4*kron(MD2N,Minv_r2);
tDelta_0=kron(I_N,MD2K)/resc^2+4*kron(MD2N,I_K/rmax^2);
itDelta=kron(I_N,derCheb_mat(2*K-1)^2)/iresc^2+4*sparse(my_kron(MD2N,convomat_coscos(iinv_r2,2*K-1)));%Properly extended to remove truncation errors later on
itDelta_0=kron(I_N,derCheb_mat(2*K-1)^2)/iresc^2+4*kron(MD2N,speye(2*K-1))/irmax^2;
iD2Psi=kron(MD2N,speye(2*K-1));

%Matrices of the scalar product 
Mat_scal=scalFourCheb_mat(K,N);
iMat_scal=iscalFourCheb_mat(2*K-1,N);

%Dirichlet boundary condition in r
Mat_DBC=[zeros(2,K-2);speye(K-2)];
Mat_DBC(1,1:2:end)=-2;
Mat_DBC(2,2:2:end)=-1;
Mat_DBC_tens=kron(I_N,Mat_DBC);
iMat_DBC_tens=intval(Mat_DBC_tens);

Trunc=[speye(K-2),zeros(K-2,2)];
Trunc_tens=kron(I_N,Trunc);

Reset_for_BC=[zeros(K-2,2),speye(K-2)];
Reset_for_BC_tens=kron(I_N,Reset_for_BC);

Pad=kron(I_N,[I_K;zeros(K-1,K)]);


%% Structure for tD
operator_tD.type='tD';

operator_tD.I_K=I_K;
operator_tD.Id=I_tens;
operator_tD.MD2K=MD2K;

operator_tD.H=-tDelta;
operator_tD.H0=-tDelta_0;
operator_tD.Mat_scal=Mat_scal;
operator_tD.Mat_BC=Mat_DBC;
operator_tD.Mat_BC_tens=Mat_DBC_tens;
operator_tD.Mat_BC_0_tens=Mat_DBC_tens;
operator_tD.Reset_for_BC_tens=Reset_for_BC_tens;
operator_tD.Trunc=Trunc;
operator_tD.Trunc_tens=Trunc_tens;
operator_tD.Pad=Pad;
operator_tD.coef_C0_tens=repmat([1 2*ones(1,2*(K-1))],[1,2*N-1]);

operator_tD.iH=-itDelta;
operator_tD.iH0=-itDelta_0;
operator_tD.iMat_scal=iMat_scal;
operator_tD.iMat_BC_tens=iMat_DBC_tens;
operator_tD.iMat_BC_0_tens=iMat_DBC_tens;

operator_tD.inv_r2_err=inv_r2_err;
operator_tD.icoef_inv_r2=4*iD2Psi;

%% All the matrix representations of the operators involved in L and L^*
%Some of them were already defined, but have to be "extended" because there
%are nonlinear terms of higher degree in L than in tD, and therefore we
%need more modes to prevent truncation errors.

I_K_ext=speye(3*K-2);
I_N_ext=speye(2*N+3);

%vectors 1 and r in Chebyshev (on [-1,1])
ie0K=intval(zeros(3*K-2,1));
ie0K(1)=1;
iM0K=sparse(convomat_coscos(ie0K));

%derivatives in psi
MDN_ext=1i*spdiags((-N-1:N+1)',0,2*N+3,2*N+3);
MDN_tens_ext=kron(MDN_ext,I_K_ext);
MD2N_ext=-spdiags((-N-1:N+1)'.^2,0,2*N+3,2*N+3);

%derivatives in r
MDK_ext=derCheb_mat(3*K-2);
MDK_tens_ext=kron(I_N_ext,MDK_ext);
MD2K_ext=MDK_ext^2;

%the operator \tilde{\Delta} 
itDelta=kron(I_N_ext,MD2K_ext)/iresc^2+4*sparse(my_kron(MD2N_ext,convomat_coscos(iinv_r2,3*K-2)));%Properly extended to remove truncation errors later on

%f% 
f=alpha*r-a*r3+sigma^2/2*inv_r;
Mf=kron(I_N,convomat_coscos(f));
norm_f2=(abs(f(1))+2*sum(abs(f(2:end))))^2;
ir3=iMr*ir2;
i_f=ialpha*ir-ia*ir3+isigma^2/2*iinv_r;
iMf=my_kron(I_N_ext,convomat_coscos(i_f,3*K-2));
inorm_f2=(abs(i_f(1))+2*sum(abs(i_f(2:end))))^2;
% figure
% plotCheb(f,rmin,rmax);


%g
e0N=zeros(2*N-1,1);
e0N(N)=1;
cos1N=zeros(2*N-1,1);
cos1N(N-1)=1/2;
cos1N(N+1)=1/2;
ie0N=zeros(2*N+3,1);
ie0N(N+2)=1;
icos1N=zeros(2*N+3,1);
icos1N(N+1)=1/2;
icos1N(N+3)=1/2;
Mg=2*kron(convomat_exp(b*e0N+sqrt(a^2+b^2)*cos1N),Mr2);
norm_r2g2=(rmax^3*(b+sqrt(a^2+b^2)))^2;
iMg=2*my_kron(sparse(convomat_exp(ib*ie0N+sqrt(ia^2+ib^2)*icos1N)),convomat_coscos(ir2,3*K-2));
inorm_r2g2=(irmax^3*(ib+sqrt(ia^2+ib^2)))^2;
% figure
% plotFourCheb(g,rmin,rmax);

%h=Drf+DPsig
isin1N=zeros(2*N+3,1);
isin1N(N+1)=1i/2;
isin1N(N+3)=-1i/2;
isin1N=sparse(isin1N);
iMh=my_kron(I_N_ext,convomat_coscos(ialpha*e0K-3*ia*ir2-isigma^2/2*iinv_r2,3*K-2))...
    -2*my_kron(convomat_exp(sqrt(ia^2+ib^2)*isin1N),convomat_coscos(ir2,3*K-2));

%h0=inf(df/dr+dg/dpsi)=inf h
C1=3*a+2*sqrt(a^2+b^2);
C2=sigma^2/2;
h0=alpha-max(C1*rmax^2+C2/rmax^2,C1*rmin^2+C2/rmin^2);
iC1=3*ia+2*sqrt(ia^2+ib^2);
iC2=isigma^2/2;
ih0=inf(ialpha-max(iC1*irmax^2+iC2/irmax^2,iC1*irmin^2+iC2/irmin^2));

%h1 = sup h, ||h||_\infty = max(|h0|,|h1|)
iC1=3*ia-2*sqrt(ia^2+ib^2);
iC2=isigma^2/2;
if iC1.sup<0 %C1<0
    ih1=sup(ialpha-iC1*irmax^2-iC2/irmax^2);
elseif iC1.inf>0 %C1>0
    ir_extr=(iC2/iC1)^(1/4);
    if (ir_extr.sup >= irmin) && (ir_extr.inf <= irmax)
        ih1=sup(ialpha-iC1*ir_extr^2-iC2/ir_extr^2);
    else
        ih1=sup(ialpha-min(iC1*irmax^2+iC2/irmax^2,iC1*irmin^2+iC2/irmin^2));
    end
else %C1 may be equal to 0
    ih1=sup(ialpha-iC1*infsup(irmin.inf,irmax.sup)^2-iC2/irmax^2);
end
inorm_h2=intval(max(abs(ih0),abs(ih1)))^2;
norm_h2=inorm_h2.mid;


%V, Vstar
V=Mf/resc*MDK_tens+Mg*MDN_tens;
Vstar=-(MDK_tens*Mf/resc+MDN_tens*Mg);
iV=iMf/iresc*MDK_tens_ext+iMg*MDN_tens_ext;
iVstar=-(MDK_tens_ext*iMf/iresc+MDN_tens_ext*iMg);


%L, Lstar, cL
L=sigma^2/2*tDelta+V;
Lstar=sigma^2/2*tDelta+Vstar;
iL=isigma^2/2*itDelta+iV;
iLstar=isigma^2/2*itDelta+iVstar;

%% Structure for L
operator_L.L=L;
operator_L.Lstar=Lstar;
operator_L.Mat_BC_tens=Mat_DBC_tens;
operator_L.Trunc_tens=Trunc_tens;
operator_L.Reset_for_BC_tens=Reset_for_BC_tens;

operator_L.norm_f2=norm_f2;
operator_L.norm_r2g2=norm_r2g2;
operator_L.norm_h2=norm_h2;
operator_L.h0=h0;

operator_L.itDelta=itDelta;
operator_L.iL=iL;
operator_L.iLstar=iLstar;
operator_L.iMat_BC_tens=iMat_DBC_tens;

operator_L.iMf=iMf;
operator_L.iMg=iMg;
operator_L.iMh=iMh;

operator_L.inorm_f2=inorm_f2;
operator_L.inorm_r2g2=inorm_r2g2;
operator_L.inorm_h2=inorm_h2;
operator_L.ih0=ih0;
operator_L.ih1=ih1;

operator_L.I_N_ext=I_N_ext;
operator_L.iM0K=iM0K;
operator_L.isin1N=isin1N;
operator_L.ir2=ir2;
