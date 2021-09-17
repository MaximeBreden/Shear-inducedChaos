function [eigenval_s,eigenvect_s,ilower_eig_s]=selecting_eigs(eigenval_tD,eigenvect_tD,ieigenval_tD,op,eig_target)

%Here we select a close to minimal number of eigenvalue to start the
%homotopy (starting with extra eigenvalues is fine theoretically, but
%computationally inefficient).

M=length(eigenval_tD);
eigenval_S0=[op.s2*eigenval_tD.^2-op.s1*eigenval_tD+op.s0;op.slambda];
eigenvect_S0=[[eigenvect_tD;zeros(1,M)],[zeros(size(eigenvect_tD,1),1);1]];
M=M+1;%we just added the eigenvalue s_\lambda and the corresponding eigenvector [0,...,0,1]
[eigenval_S0,ind]=sort(eigenval_S0);
eigenvect_S0=eigenvect_S0(:,ind);
if max(eigenval_S0)<eig_target
    error('We do not have enough eigenvalues of tDelta, the second homotopy will not be able to reach s=1.')
end
m=M;
while m>1 && ( eigenval_S0(m-1)>max(eig_target+100,1.1*eig_target) || abs((eigenval_S0(m)-eigenval_S0(m-1))/eigenval_S0(m-1))<10^-8 ) %removing eigenvalues if we have too much margin, or if they are clustered
    m=m-1;
end

ieigenval_S0=[op.is2*ieigenval_tD.^2-op.is1*ieigenval_tD+op.is0;op.islambda];
ieigenval_S0=ieigenval_S0(ind);

eigenval_s=eigenval_S0(1:m);
eigenvect_s=eigenvect_S0(:,1:m);
ilower_eig_s=sort(ieigenval_S0.inf,'ascend');
ilower_eig_s=ilower_eig_s(1:m);

%Sanity check
eig_min=op.s0-op.s1^2/(2*op.s2);
[~,eigenval_S0_bis]=eigs(op.Trunc_tens*(op.H0-eig_min*op.Id)*op.Mat_BC_0_tens,op.Trunc_tens*op.Mat_BC_0_tens,m,'SM');
eigenval_S0_bis=sort(real(diag(eigenval_S0_bis))+eig_min,'ascend');
diff=max(abs((eigenval_S0(1:m)-eigenval_S0_bis))./abs(eigenval_S0(1:m)));
if diff>10^-4
    warning('Something may be wrong with the computation of the eigenvalues of S0')
    [eigenval_S0(1:m) eigenval_S0_bis abs((eigenval_S0(1:m)-eigenval_S0_bis))./abs(eigenval_S0(1:m))]
    pause
end

disp('Eigenvalues at the start of the homotopy:')
eigenval_s