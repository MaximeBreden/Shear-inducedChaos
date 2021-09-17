function [eigenval,eigenvect,s]=homotopy_breakpoint_check(s,delta_s,nu_min,nu_max,tol_nu,M,eigenval,eigenvect,op,delta_s_min,show)

%Here we find a breakpoint s_{k+1} for the homotopy, starting from the 
%initial guess s+delta_s and using a bissection algorithm. 

eig_max=max(eigenval);
if not(eig_max<nu_max)
    warning('We are already too far up')
    return
end

close_enough=(abs(nu_max-eig_max)/abs(nu_max)<tol_nu);
threshold_reached=0;
while eig_max<nu_max && s<1 && not(close_enough) && delta_s>delta_s_min
    s_old=s;
    eigenvect_old=eigenvect;
    eigenval_old=eigenval;
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
    if eig_max>=nu_max        
        threshold_reached=1;
        s=s_old;
        eigenvect=eigenvect_old;
        eigenval=eigenval_old;
        eig_max=max(eigenval);
        delta_s=delta_s/2;
    elseif threshold_reached
        delta_s=delta_s/2;
    end
    close_enough=(abs(nu_max-eig_max)/abs(nu_max)<tol_nu);
end

if eig_max>=nu_max
    s=s_old;
    eigenvect=eigenvect_old;
    eigenval=eigenval_old;
end

