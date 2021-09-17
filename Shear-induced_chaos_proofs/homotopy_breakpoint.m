function [eigenval,eigenvect,s]=homotopy_breakpoint(s,delta_s,nu_min,nu_max,tol_nu,M,eigenval,eigenvect,op,threshold_reached,show)

%Here we find a breakpoint s_{k+1} for the homotopy, starting from the 
%initial guess s+delta_s and using a bissection algorithm. 

eig_max=max(eigenval);
test_nu=eig_max<nu_max;
if not(test_nu)
    disp('We are already too far up')
    return
end

close_enough=(abs(nu_max-eig_max)/abs(nu_max)<tol_nu);
while test_nu && s<1 && (not(threshold_reached) || (threshold_reached && not(close_enough)))
    s_old=s;
    eigenvect_old=eigenvect;
    eigenval_old=eigenval;
    eig_max_old=eig_max;
    s=min(s+delta_s,1);    
    if show
        fprintf('\ns=%f\n',s)
    end
    Mat_BC_s_tens=s*op.Mat_BC_tens+(1-s)*op.Mat_BC_0_tens;
    A=op.Trunc_tens*(s*op.H+(1-s)*op.H0-nu_min*op.Id)*Mat_BC_s_tens;
    B=op.Trunc_tens*Mat_BC_s_tens;
    [eigenvect,eigenval,flag]=eigs(A,B,M,'SM');
    if flag~=0
        fprintf('eigs may not have converged')
    end
    eigenval=nu_min+real(diag(eigenval));
    eig_max=max(eigenval);
    test_nu=eig_max<nu_max;
    close_enough=(abs(nu_max-eig_max)/abs(nu_max)<tol_nu);
    if threshold_reached && test_nu
        delta_s=delta_s/2;
    end
end

if not(test_nu)
    if abs(nu_max-eig_max_old)/abs(nu_max)>tol_nu
        threshold_reached=1;
        [eigenval,eigenvect,s]=homotopy_breakpoint(s_old,delta_s/2,nu_min,nu_max,tol_nu,M,eigenval_old,eigenvect_old,op,threshold_reached,show);
     else
        s=s_old;
        eigenvect=eigenvect_old;
        eigenval=eigenval_old;
    end
end


