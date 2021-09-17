%%%% Non rigorous computation of the Lyapunov exponent %%%%
%%%% Not to be fully trusted without validation, especially for sigma small %%%%

clear variables
close all
clc
warning('on')

%Parameters of the original SDE (alpha, a, b, sigma)
alpha=1;
a=1;
b=3.5;

SDE.alpha=alpha;
SDE.a=a;
SDE.b=b;

%Domain (rmin, rmax)
rmin=0.75;
rmax=1.25;
resc=(rmax-rmin)/2;
Domain.rmin=rmin;
Domain.rmax=rmax;
Domain.resc=resc;


K=30;%Chebyshev in r
N=30;%Fourier in psi

% [Lyap,eta,phi]=approx_Lyap(SDE,Domain,K,N);
% sigma
% Lyap
% figure
% plotFourCheb(eta,rmin,rmax);
% title(['\eta, for \sigma =',num2str(sigma),', b =',num2str(b)])
% figure
% plotFourCheb(phi,rmin,rmax);
% title(['\phi, for \sigma =',num2str(sigma),', b =',num2str(b)])
% drawnow


Tab_sigma=0.25:0.1:1.75;
Tab_Lyap=0*Tab_sigma;
for l=1:length(Tab_sigma)
    sigma=Tab_sigma(l)
    SDE.sigma=sigma;
    [Lyap,~,~]=approx_Lyap(SDE,Domain,K,N);
    Tab_Lyap(l)=Lyap;
end

figure
plot(Tab_sigma,Tab_Lyap,'*k')
xlabel('\sigma')
ylabel('\Lambda_c')
title(['\alpha =',num2str(alpha),', a =',num2str(a),', b =',num2str(b),...
       ', r_{min} =',num2str(rmin),', r_{max} =',num2str(rmax)])
   

   
