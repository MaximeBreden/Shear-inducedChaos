function [delta_s,correc]=delta_s_predictor(s,nu_max,tol_nu,op,small_nb_eigs,pert,delta_s_min,correc,show)

%Cheap approximate computation/guess of delta_s = s_{k+1}-s_{k}. We shift 
%the spectrum so that the breakpoint corresponds to one eigenvalue crossing
%0, and only look at a few eigenvalues around 0. Of course since we work 
%with a finite stepsize for s, we might miss a crossing if some eigenvalues
%change much more rapidly than others with respect to s.

shift=nu_max;
Mat_BC_s_tens=s*op.Mat_BC_tens+(1-s)*op.Mat_BC_0_tens;
A=op.Trunc_tens*(s*op.H+(1-s)*op.H0-shift*op.Id)*Mat_BC_s_tens;
B=op.Trunc_tens*Mat_BC_s_tens;
[~,eigenval,flag]=eigs(A,B,small_nb_eigs,'SM');
if flag~=0
    fprintf('eigs may not have converged')
end
eigenval=real(diag(eigenval));
eigenval=sort(eigenval,'ascend');
nb_pos=sum(eigenval>0);
nb_neg=small_nb_eigs-nb_pos;

if nb_neg<2
    delta_s=delta_s_predictor(s,nu_max,tol_nu,op,small_nb_eigs+2,pert,delta_s_min,correc,show);
    return
end
eig_max=eigenval(nb_neg);
if eig_max>=nu_max
    warning('We are already too far up')
    delta_s=delta_s_min;
    return
end

s_pert=s+pert;
Mat_BC_s_tens=s_pert*op.Mat_BC_tens+(1-s_pert)*op.Mat_BC_0_tens;
A=op.Trunc_tens*(s_pert*op.H+(1-s_pert)*op.H0-shift*op.Id)*Mat_BC_s_tens;
B=op.Trunc_tens*Mat_BC_s_tens;
[~,eigenval_pert,flag]=eigs(A,B,small_nb_eigs,'SM');
if flag~=0
    fprintf('eigs may not have converged')
end
eigenval_pert=real(diag(eigenval_pert));
eigenval_pert=sort(eigenval_pert,'ascend');

%If eigs does not really give the smallest eigenvalues, which apparently
%happens sometimes
diff1=(eigenval_pert-eigenval);
diff2=(eigenval_pert(1:end-1)-eigenval(2:end));
diff3=(eigenval_pert(2:end)-eigenval(1:end-1));
n1=norm(diff1,1);
n2=norm(diff2,1);
n3=norm(diff3,1);
if n1<=min(n2,n3)
    dlambda_ds=(diff1)/pert;
    ds=max((1-correc)*min(-eigenval(1:nb_neg)./dlambda_ds(1:nb_neg)),100*delta_s_min);
else
    if n2<=n3
        dlambda_ds=(diff2)/pert;
        ds=max((1-correc)*min(-eigenval(2:nb_neg)./dlambda_ds(1:nb_neg-1)),100*delta_s_min);
    else
        dlambda_ds=(diff3)/pert;
        ds=max((1-correc)*min(-eigenval(1:nb_neg)./dlambda_ds(1:nb_neg)),100*delta_s_min);
    end
end
ds=min(ds,(1-s)/10);

sf=1-10^-8;%security factor, to account for some small numerical errors
close_enough=abs(eig_max/nu_max)<tol_nu;
s0=s;
threshold_reached=0;
it=0;
while eig_max<0 && s<1 && not(close_enough) && ds>=delta_s_min
    if it==1 && ds>delta_s_min
        if threshold_reached
            correc=min(5*correc,10^-1);
        else
            correc=max(correc/5,10^-6);
        end
    end
    s_old=s;
    eigenval_old=eigenval;
    nb_neg_old=nb_neg;
    nb_pos_old=nb_pos;
    eig_max_old=eig_max;
    s=min(s+ds,1);    
    if show
        fprintf('\ns=%f\n',s)
    end
    Mat_BC_s_tens=s*op.Mat_BC_tens+(1-s)*op.Mat_BC_0_tens;
    A=op.Trunc_tens*(s*op.H+(1-s)*op.H0-shift*op.Id)*Mat_BC_s_tens;
    B=op.Trunc_tens*Mat_BC_s_tens;
    [~,eigenval,flag]=eigs(A,B,small_nb_eigs,'SM');
    if flag~=0
        fprintf('eigs may not have converged')
    end
    eigenval=real(diag(eigenval));
    eigenval=sort(eigenval,'ascend');
    nb_pos=sum(eigenval>=0);
    nb_neg=small_nb_eigs-nb_pos;
    if nb_neg==0
        threshold_reached=1;
        s=s_old;
        eigenval=eigenval_old;
        nb_neg=nb_neg_old;
        nb_pos=nb_pos_old;
        eig_max=eig_max_old;
        ds=ds/2;
    else
        eig_max=eigenval(nb_neg);
        if nb_pos>0
            if nb_pos>nb_pos_old || (eig_max<sf*eig_max_old && s-s_old>eps) || max(eigenval(nb_neg+1:end)<sf*eigenval_old(nb_neg_old+1:nb_neg_old+nb_pos))
                threshold_reached=1;
                s=s_old;
                eigenval=eigenval_old;
                nb_neg=nb_neg_old;
                nb_pos=nb_pos_old;
                eig_max=eig_max_old;
                ds=ds/2;
            elseif threshold_reached
                ds=ds/2;
            end
        else
            if eig_max<sf*eig_max_old 
                threshold_reached=1;
                s=s_old;
                eigenval=eigenval_old;
                nb_neg=nb_neg_old;
                nb_pos=nb_pos_old;
                eig_max=eig_max_old;
                ds=ds/2;
            elseif threshold_reached
                ds=ds/2;
            end
        end
    end
    close_enough=abs(eig_max/nu_max)<tol_nu;
    it=it+1;
end

delta_s=s-s0;

