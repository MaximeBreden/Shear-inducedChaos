%%%% rigorous computation of the Lyapunov exponent %%%%

% You need to start Intlab before running this code 

clear variables
close all
clc
warning('on')
intvalinit('FastIVmult')

tic;

%Parameters of the original SDE (alpha, a, b, sigma)
ialpha=intval('1');
ia=intval('1');
ib=intval('3.5');
isigma=intval('1.2');

SDE.ialpha=ialpha;
SDE.ia=ia;
SDE.ib=ib;
SDE.isigma=isigma;
SDE.alpha=mid(SDE.ialpha);
SDE.a=mid(SDE.ia);
SDE.b=mid(SDE.ib);
SDE.sigma=mid(SDE.isigma);

%Domain (rmin, rmax)
irmin=intval('0.75');
irmax=intval('1.25');
rmin=irmin.mid;
rmax=irmax.mid;
resc=(rmax-rmin)/2;
iresc=(irmax-irmin)/2;

Domain.irmin=irmin;
Domain.irmax=irmax;
Domain.rmin=rmin;
Domain.rmax=rmax;
Domain.resc=resc;
Domain.iresc=iresc;

%A couple of preliminary computations, including a rough estimation of the largest eigenvalue of \tilde{\Delta} we will need to rigorously enclose
[lb1,lb2,eta_L1,eta_L2]=rough_lb_nbeig_tD(SDE,Domain);
drawnow
lb=max(lb1,lb2);
lb_margin=max(lb+10,1.1*lb);%to start will a little bit of margin (can be changed if needed)

%Computational parameters
K=30;%Chebyshev in r
N=30;%Fourier in psi

%Weights for the norms (since we always ended up using xi0 = xi1 = 1, these
%weights have been removed from the paper)
ixi0=intval('1');
xi0=ixi0.mid;
ixi1=intval('1');
xi1=ixi1.mid;

weights.ixi0=ixi0;
weights.ixi1=ixi1;
weights.xi0=xi0;
weights.xi1=xi1;
weights.eta_L1=eta_L1;
weights.eta_L2=eta_L2;

%for saving the results
str=['proof_[',num2str(rmin),',',num2str(rmax),']_b=',num2str(SDE.b),'_sigma=',num2str(SDE.sigma),'.mat'];

%% Initialization of all the matrices needed for the first homopty
fprintf('\n\nInitialization of some operators related to tD and L ...\n')
[op_tD,op_L,approx_r]=initialize_tD_L(SDE,Domain,K,N);

%% Finding the number M of eigenvalues of the base problem \tilde{\Delta}^{(0)} that we need
fprintf('\n\nFinding M...\n')
k_max=ceil((rmax-rmin)/pi*sqrt(lb_margin))+1;
n_max=ceil(2/rmax*sqrt(lb_margin))+1;
[k,n]=meshgrid(1:k_max,-n_max:n_max);

ipi=intval('pi');
iphi=(ipi/(irmax-irmin))^2*k.^2+(2/irmax)^2*n.^2;

phi=(pi/(rmax-rmin))^2*k.^2+(2/rmax)^2*n.^2;
M=sum(sum(phi<=lb_margin));
[~,ind_phi]=sort(reshape(phi,[numel(phi),1]));
ind_phi_remainder=ind_phi(M+1:end);
ind_phi=ind_phi(1:M);

if max(abs(n(ind_phi)))>=N
    fprintf('\nNot enough Fourier mode, increase N to at least %d\n',max(abs(n(ind_phi)))+1)
    return
end

%Taking care of the fact that the "last" eigenvalue is potentially of multiplicity 2
if k(ind_phi(end))==k(ind_phi_remainder(1)) && n(ind_phi(end))==-n(ind_phi_remainder(1))%If the M-th eigenvalue is double and only one was kept
    ind_phi_remainder=ind_phi_remainder(2:end);
elseif k(ind_phi(end))==k(ind_phi(end-1)) && n(ind_phi(end))==-n(ind_phi(end-1))%If the M-th eigenvalue is double and both were kept
    ind_phi=ind_phi(1:end-1);
    M=M-1;
end

%% Computation of the M smallest eigenvalues, and associated eigenvectors, of the base problem \tilde{\Delta}^{(0)} = -(d_rr  + 4/rmax^2 d_psipsi )
eigenval_0=phi(ind_phi);%The M smallest eigenvalues
ieigenval_0=iphi(ind_phi);

%Let's rigorously check that we indeed have all the smallest eigenvalues
smallest_eig_not_kept=min([min(iphi.inf(ind_phi_remainder)),(ipi/(irmax-irmin))^2*(k_max+1)^2,(2/irmax)^2*(n_max+1)^2]);
if not(max(ieigenval_0.sup)<smallest_eig_not_kept)
    error('Not able to guarantee the we indeed start with the M smallest eigenvalues')
end

disp('The starting eigenvalues are:')
disp(eigenval_0)
disp(['M = ',num2str(M)])

%Constructing the associated eigenvectors (does not have to be rigorous),
k=k(ind_phi);
n=n(ind_phi);
I_K=speye(K);
eigenvect_k=zeros(K,max(k));%eigenvectors in the r variable only
for l=1:max(k)
    lambda_l=(pi/(rmax-rmin))^2*l^2;
    eigenvect_l=op_tD.Mat_BC*null(op_tD.Trunc*(-full(op_tD.MD2K)/resc^2-lambda_l*op_tD.I_K)*op_tD.Mat_BC);
    if isempty(eigenvect_l)
        [eigenvect_l,~]=eigs(op_tD.Trunc*(-op_tD.MD2K/resc^2-lambda_l*op_tD.I_K)*op_tD.Mat_BC,1,'SM');
        eigenvect_l=op_tD.Mat_BC*eigenvect_l;
    end        
    eigenvect_k(:,l)=eigenvect_l;
end
eigenvect_n=zeros(2*N-1,max(n)-min(n)+1);%eigenvectors in the theta variable only
for m=1:M
    if n(m)==0
        eigenvect_n(N,n(m)-min(n)+1)=1;
    elseif n(m)>0
        eigenvect_n(N+[-n(m) n(m)],n(m)-min(n)+1)=1/2;
    else
        eigenvect_n(N+[-n(m) n(m)],n(m)-min(n)+1)=[-1i/2 1i/2];
    end
end

eigenvect_0=zeros(K*(2*N-1),M);%the eigenvectors associated to the M smallest eigenvalues
for m=1:M
    eigenvect_0(:,m)=reshape(kron(transpose(eigenvect_n(:,n(m)-min(n)+1)),eigenvect_k(:,k(m))),[K*(2*N-1),1]);
end

%Sanity check
[~,eigenval_00]=eigs(op_tD.Trunc_tens*op_tD.H0*op_tD.Mat_BC_tens,op_tD.Trunc_tens*op_tD.Mat_BC_tens,M,'SM');
eigenval_00=real(diag(eigenval_00));
if max(abs(eigenval_0-eigenval_00)./abs(eigenval_0))>10^-3
    [eigenval_0,eigenval_00,abs(eigenval_0-eigenval_00)./abs(eigenval_0)]
    warning('You may not have enough modes to represent the eigenvectors accurately')
    %pause
end

%% First homotopy to get the eigenvalues of tD 
fprintf('\n\nFirst homotopy to get the eigenvalues of tD\n')
s=0;
steps_s=0;
margin_tD=margin_for_homotopy(eigenval_0);%This parameter controls how the breakpoints s_k for the homotopy are chosen 
fprintf('\nmargin_tD = %f\n',margin_tD)
tol_simple_eig=10^-7;%For the validation of simple VS clustered eigenvalues
show_s=0;%Put to 1 to display what the value of s is during the whole process
show_enclosure=0;%Put to 1 to display the obtained enclosures at each homotopy breakpoint

%The homotopy method (algorithm 1 in the paper)
[ilower_eig_tD,iupper_eig_tD,eigenvect_tD,steps_s]=homotopy(s,steps_s,op_tD,approx_r,eigenval_0,eigenvect_0,ieigenval_0.inf,margin_tD,tol_simple_eig,show_s,show_enclosure);

ieigenval_tD=infsup(ilower_eig_tD,iupper_eig_tD);%rigorous enclosure of the first eigenvalues of tD
eigenval_tD=ieigenval_tD.mid;

fprintf('\nlb = %g, largest validated eigenvalue of tD = %g\n',lb,max(eigenval_tD))

%Rough estimate of whether we have enough eigenvalues of tD
if lb>max(eigenval_tD)
    warning('We may not have enough eigenvalues of tD')
end

str_tD=['logtD_',str];
time=toc;
save(str_tD,'SDE','Domain','K','N','steps_s','ieigenval_tD','time')

%% Getting the principal eigenvalue/eigenvector of L and L*
fprintf('\n\nComputing the principal eigenvalue/eigenvector of L and L^*\n')

[PEVect_L,PEVal_L]=eigs(op_L.Trunc_tens*op_L.L*op_L.Mat_BC_tens,op_L.Trunc_tens*op_L.Mat_BC_tens,1,'SM');
PEVect_L=op_L.Mat_BC_tens*PEVect_L;
PEVect_L=reshape(PEVect_L,[K,2*N-1]);
PEVect_L=symmetrize(PEVect_L);
PEVect_L=PEVect_L/real(sqrt(scalFourCheb(PEVect_L,PEVect_L)));
PEVal_L=real(PEVal_L);

[PEVect_Lstar,PEVal_Lstar]=eigs(op_L.Trunc_tens*op_L.Lstar*op_L.Mat_BC_tens,op_L.Trunc_tens*op_L.Mat_BC_tens,1,'SM');
PEVect_Lstar=op_L.Mat_BC_tens*PEVect_Lstar;
PEVect_Lstar=reshape(PEVect_Lstar,[K,2*N-1]);
PEVect_Lstar=symmetrize(PEVect_Lstar);
PEVect_Lstar=PEVect_Lstar/real(sqrt(scalFourCheb(PEVect_Lstar,PEVect_Lstar)));
PEVal_Lstar=real(PEVal_Lstar);

%Changing the sign of the eigenvectors if they happen to be negative
%instead of positive
figure
val=plotFourCheb(PEVect_L);
if val<-10*eps
    PEVect_L=-PEVect_L;
end
close gcf
figure
val=plotFourCheb(PEVect_Lstar);
if val<-10*eps
    PEVect_Lstar=-PEVect_Lstar;
end
close gcf

%% Initialization of the remaining operators needed for the second homotpy (for L)
fprintf('\n\nInitialization of the remaining operators needed for the second homotpy (for L) ...\n')

tlambda=PEVal_L;
tU=PEVect_L;
op_S1=initialize_S1(SDE,Domain,tU,tlambda,op_tD,op_L,approx_r,weights,K,N);

%Rigorously checking that we have enough eigenvalues of tD
ieig_max=intval(max(ieigenval_tD.inf));
if not(op_S1.is2*ieig_max^2-op_S1.is1*ieig_max+op_S1.is0>0)
    error('We do not have enough eigenvalues of tD to guarantee that we have all the smallest eigenvalues of SS^*')
end

[~,eig_target]=eigs(op_S1.Trunc_tens*op_S1.H*op_S1.Mat_BC_tens,op_S1.Trunc_tens*op_S1.Mat_BC_tens,1,'SM');
eig_target=real(eig_target);
fprintf('\nThe final eigenvalue that should be the smallest and that we want to validate:\neig_target = %d\n\n',eig_target)

[eigenval_s,eigenvect_s,ilower_eig_s]=selecting_eigs(eigenval_tD,eigenvect_tD,ieigenval_tD,op_S1,eig_target);

%% Homotopy to get k0
fprintf('\n\nSecond homotopy to get k0 (for L)\n')

s=0;
steps_s=0;
margin_S1=margin_for_homotopy(eigenval_s);%This parameter controls how the breakpoints s_k for the homotopy are chosen 
fprintf('\nmargin_S1 = %f\n',margin_S1)
tol_simple_eig=10^-7;%For the validation of simple VS clustered eigenvalues
show_s=0;%Put to 1 to display what the value of s is during the whole process
show_enclosure=0;%Put to 1 to display the obtained enclosures at each homotopy breakpoint

%The homotopy method (algorithm 1 in the paper)
[ilower_eig_s,~,~,steps_s]=homotopy(s,steps_s,op_S1,approx_r,eigenval_s,eigenvect_s,ilower_eig_s,margin_S1,tol_simple_eig,show_s,show_enclosure);

ikappa0=1/intval(ilower_eig_s(1));
fprintf('\n kappa0 = %f\n',ikappa0.mid)

%% from kappa_0 to kappa
iCV=op_S1.iCV;
itu=op_S1.itu;
itlambda=op_S1.itlambda;
inorm_tu2=op_S1.inorm_tu2;

ikappa1=1/isigma^2*(iCV*ikappa0+sqrt(iCV^2*ikappa0^2+2*isigma^2*ikappa0*(1+ikappa0*max(-itlambda,0))));
fprintf('\n kappa1 = %f\n',ikappa1.sup)

ikappa2=2/isigma^2*(1+iCV*ikappa1+abs(itlambda)*ikappa0);
fprintf('\n kappa2 = %f\n',ikappa2.sup)

itheta=intval('1');

ixi2=1/16*(isigma^4/(2*inorm_tu2));

ikappa=1/sqrt(ixi0)*sqrt((ixi0*ikappa0^2+ixi1*ikappa1^2+(1+itheta)*ixi2*ikappa2^2)/(1-(1+itheta)/itheta*ixi2*(2*inorm_tu2/isigma^2)^2));
fprintf('\n kappa = %f\n',ikappa.sup)

%% Checking that we have a contraction, rigorous validation of the eigenpair
igamma=intval('1');
iMat_scal_X=op_S1.iMat_scal;

iresidue=[op_L.iL*itu-itlambda*itu;ixi0*inorm_tu2-1];
inorm_residue=sqrt(max(real(iresidue'*iMat_scal_X*iresidue),0));

iresidue_u=iresidue(1:end-1);
ierr_residue_r2=[(2*isigma^2*op_S1.iD2Psi)*iresidue_u;0];
ierr_residue_r=[(isigma^2/2*op_S1.iDr)*iresidue_u;0];
ierr_residue=approx_r.inv_r2_err*sqrt(real(ierr_residue_r2'*iMat_scal_X*ierr_residue_r2))...
             +approx_r.inv_r_err*sqrt(real(ierr_residue_r'*iMat_scal_X*ierr_residue_r));

idelta=inorm_residue+ierr_residue;
fprintf('\n delta = %e\n',idelta.sup)

idisc=1-2*ikappa^2*igamma*idelta;
fprintf('\n discriminant = %f\n',idisc.sup)

str_S1=['logS1_',str];
time=toc;
save(str_S1,'SDE','Domain','K','N','steps_s','ilower_eig_s','idelta','ikappa','time')

if idisc>0
    irho_min=(1-sqrt(idisc))/(ikappa*igamma);
    irho_max=1/(ikappa*igamma);
    fprintf('\n [rho_min,rho_max) = [%e,%e)\n',irho_min.sup,irho_max.inf)
else
    error('Unable to prove that we have a fixed point, the bounds for delta and/or kappa are not sharp enough')
end

ieta=itu;
irho_eta_min=irho_min;
irho_eta_max=irho_max;

%% Checking that we have the correct eigenpair
%Constants for the embeddings 
igamma1=intval('1.1548');
igamma2=intval('0.22361');
il1=2*ipi;
il2=irmax-irmin;
iC0=sqrt(il1*il2);
iC1=igamma1*sqrt((il1^2+il2^2)/(3*il1*il2));
iC2=igamma2*sqrt(((il1^2+il2^2)^2+4*(il1^4+il2^4)/3)/(3*il1*il2));
im1=max(irmax/2,1);
im2=max(irmax^2/4,(irmin^2+irmax^2)/(2*irmin^2));
iUpsilon_X_C0=sqrt(6*ipi*(irmax-irmin))*max([iC0/sqrt(ixi0),iC1*im1/sqrt(ixi1),iC2*im2/sqrt(ixi2)]);
iUpsilon_Xrad_C1=sqrt((rmax-rmin)/tanh(irmax-irmin))*max(1/sqrt(ixi0),1/sqrt(ixi1));

%C^0 error estimates
ieps_eta=iUpsilon_X_C0*irho_min;
ieps_der_eta=iUpsilon_Xrad_C1*irho_min;

%Rigorous estimates for the numerical eigenvector
ieta_rad=real(op_S1.itU(1:K,N+2));%1D numerical eigenvector eta
ider_eta_rad=derCheb(ieta_rad)/iresc;
ider2_eta_rad=derCheb(ider_eta_rad)/iresc;
ider2_eta_rad_max=[1 2*ones(1,K-1)]*abs(ider2_eta_rad);%An upper bound of the second derivative of eta on [rmin,rmax]
ider_eta_rad_rmin=[1 2*(-1).^(1:K-1)]*ider_eta_rad;%eta'(rmin)
if ider_eta_rad_rmin.inf<=0
    error('the derivative of eta at rmin does not have the correct sign')
end
ider_eta_rad_rmax=[1 2*ones(1,K-1)]*ider_eta_rad;%eta'(rmax)
if ider_eta_rad_rmax.sup>=0
    error('the derivative of eta at rmax does not have the correct sign')
end

ieps_rmin=(ider_eta_rad_rmin-ieps_der_eta)/ider2_eta_rad_max;
ieps_rmax=(-ider_eta_rad_rmax-ieps_der_eta)/ider2_eta_rad_max;
fprintf('\neta is at least positive on (rmin,rmin+eps_rmin) and on (rmax-eps_rmax,rmax), where eps_rmin =%f, and  eps_rmax =%f\n',ieps_rmin.inf,ieps_rmax.inf)
fprintf('\nWe still have to check the positivity on [rmin+eps_rmin,rmax-eps_rmax] = [%f,%f]\n',inf(irmin+ieps_rmin),sup(irmax-ieps_rmax))

if sup(irmax-ieps_rmax) > inf(irmin+ieps_rmin)
    nb_sub=1;
    pos=0;
    while nb_sub<10^6 && pos==0
        it=0;
        h=((irmax-ieps_rmax)-(irmin+ieps_rmin))/nb_sub;
        pts=(irmin+ieps_rmin)+(0:nb_sub)'*h;
        int=infsup(inf(pts(1:end-1)),sup(pts(2:end)))
        ieta_eval=evalCheb(ieta_rad,int,irmin,irmax)
        if min(ieta_eval>ieps_eta)
            pos=1;
        else
            nb_sub=nb_sub*10;
        end
    end
else
    pos=1;
end

if pos
    fprintf('\nPositivity rigorously checked\n')
else
    error('Unable to prove that eta is positive')
end

fprintf('\nTHE EIGENPAIR FOR L HAS BEEN RIGOROUSLY VALIDATED!\n\n')


%% Initialization of the remaining operators needed for the second homotpy (for L^*)
fprintf('\n\nInitialization of the remaining operators needed for the second homotpy (for L^*) ...\n')
tlambda=PEVal_Lstar;
tU=PEVect_Lstar;
op_S2=initialize_S2(SDE,Domain,tU,tlambda,op_tD,op_L,op_S1,approx_r,weights,K,N);

%Rigorously checking that we have enough eigenvalues of tD
ieig_max=intval(max(ieigenval_tD.inf));
if not(op_S2.is2*ieig_max^2-op_S2.is1*ieig_max+op_S2.is0>0)
    error('We do not have enough eigenvalues of tD to guarantee that we have all the smallest eigenvalues of SS^*')
end

[~,eig_target]=eigs(op_S2.Trunc_tens*op_S2.H*op_S2.Mat_BC_tens,op_S2.Trunc_tens*op_S2.Mat_BC_tens,1,'SM');
eig_target=real(eig_target);
fprintf('\nThe final eigenvalue that should be the smallest and that we want to validate:\neig_target = %d\n\n',eig_target)

[eigenval_s,eigenvect_s,ilower_eig_s]=selecting_eigs(eigenval_tD,eigenvect_tD,ieigenval_tD,op_S2,eig_target);


%% Homotopy to get k0
fprintf('\n\nSecond homotopy to get k0 (for L^*)\n')

s=0;
steps_s=0;
margin_S2=margin_for_homotopy(eigenval_s);%This parameter controls how the breakpoints s_k for the homotopy are chosen 
fprintf('\nmargin_S2 = %f\n',margin_S2)
tol_simple_eig=10^-7;%For the validation of simple VS clustered eigenvalues
show_s=0;%Put to 1 to display what the value of s is during the whole process
show_enclosure=0;%Put to 1 to display the obtained enclosures at each homotopy breakpoint

%The homotopy method (algorithm 1 in the paper)
[ilower_eig_s,~,~,steps_s]=homotopy(s,steps_s,op_S2,approx_r,eigenval_s,eigenvect_s,ilower_eig_s,margin_S2,tol_simple_eig,show_s,show_enclosure);

ikappa0=1/intval(ilower_eig_s(1));
fprintf('\n kappa0 = %f\n',ikappa0.mid)


%% from kappa_0 to kappa
iCV=op_S2.iCV;
itu=op_S2.itu;
itlambda=op_S2.itlambda;
inorm_tu2=op_S2.inorm_tu2;
inorm_hplambda=op_S2.inorm_hplambda;
ih0=op_L.ih0;

ikappa1=1/isigma^2*(iCV*ikappa0+sqrt(iCV^2*ikappa0^2+2*isigma^2*ikappa0*(1+ikappa0*max(ih0-itlambda,0))));
fprintf('\n kappa1 = %f\n',ikappa1.sup)

ikappa2=2/isigma^2*(1+iCV*ikappa1+inorm_hplambda*ikappa0);
fprintf('\n kappa2 = %f\n',ikappa2.sup)

itheta=intval('1');

ixi2=1/16*(isigma^4/(2*inorm_tu2));

ikappa=1/sqrt(ixi0)*sqrt((ixi0*ikappa0^2+ixi1*ikappa1^2+(1+itheta)*ixi2*ikappa2^2)/(1-(1+itheta)/itheta*ixi2*(2*inorm_tu2/isigma^2)^2));
fprintf('\n kappa = %f\n',ikappa.sup)

%% Checking that we have a contraction, rigorous validation of the eigenpair
igamma=intval('1');

iresidue=[op_L.iLstar*itu-itlambda*itu;ixi0*inorm_tu2-1];
inorm_residue=sqrt(max(real(iresidue'*iMat_scal_X*iresidue),0));

iresidue_u=iresidue(1:end-1);
ierr_residue_r2=[(2*isigma^2*op_S2.iD2Psi-isigma^2/2*op_S2.I_tens_ext)*iresidue_u;0];
ierr_residue_r=[(-isigma^2/2*op_S2.iDr)*iresidue_u;0];
ierr_residue=approx_r.inv_r2_err*sqrt(real(ierr_residue_r2'*iMat_scal_X*ierr_residue_r2))...
             +approx_r.inv_r_err*sqrt(real(ierr_residue_r'*iMat_scal_X*ierr_residue_r));

idelta=inorm_residue+ierr_residue;
fprintf('\n delta = %e\n',idelta.sup)

idisc=1-2*ikappa^2*igamma*idelta;
fprintf('\n discriminant = %f\n',idisc.sup)

str_S2=['logS2_',str];
time=toc;
save(str_S2,'SDE','Domain','K','N','steps_s','ilower_eig_s','idelta','ikappa','time')

if idisc>0
    irho_min=(1-sqrt(idisc))/(ikappa*igamma);
    irho_max=1/(ikappa*igamma);
    fprintf('\n [rho_min,rho_max) = [%e,%e)\n',irho_min.sup,irho_max.inf)
else
    error('Unable to prove that we have a fixed point, the bounds for delta and/or kappa are not sharp enough')
end

iphi=itu;
irho_phi_min=irho_min;
irho_phi_max=irho_max;
%% Checking that we have the correct eigenpair

pos=real(iphi'*iscalFourCheb_mat(3*K-2,N+2)*ieta)>0;

if not(pos)
    error('Unable to prove that <eta,phi> is positive')
end

fprintf('\nTHE EIGENPAIR FOR L^* HAS BEEN RIGOROUSLY VALIDATED!\n\n')
    

%% Rigorously computing the conditioned Lyapunov exponent
fprintf('\nRigorous computation of the conditioned Lyapunov exponent:\n')

iMe=my_kron(op_L.I_N_ext,ialpha*op_L.iM0K-2*ia*convomat_coscos(op_L.ir2,3*K-2))+sqrt(ia^2+ib^2)*my_kron(convomat_exp(op_L.isin1N),convomat_coscos(op_L.ir2,3*K-2));
iLyapunov=real(iphi'*iscalFourCheb_mat(3*K-2,N+2)*iMe*ieta);

coef_C0_tens=repmat([1 2*ones(1,3*(K-1))],[1,2*N+3]);
inorm_eta_C0=coef_C0_tens*abs(ieta);
inorm_phi_C0=coef_C0_tens*abs(iphi);
ieps_phi=iUpsilon_X_C0*irho_min;
inorm_e_C0=abs(ialpha)+irmax^2*(sqrt(ia^2+ib^2)+2*abs(ia));
iLyapunov_err=inorm_e_C0*(inorm_eta_C0*ieps_phi+inorm_phi_C0*ieps_eta);

iLyapunov=iLyapunov+midrad(0,iLyapunov_err.sup);
infsup(iLyapunov)

%% Saving the result
time=toc;
save(str,'SDE','Domain','K','N','ieta','iphi','irho_eta_min','irho_phi_min','iLyapunov','time')