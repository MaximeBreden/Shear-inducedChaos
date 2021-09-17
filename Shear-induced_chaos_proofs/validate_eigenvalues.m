function ieigs=validate_eigenvalues(eigs,vects,iA,iB,tol_simple_eig,all_eigs)

%Computes a rigorous enclosure of each eigenvalue of the generalized
%eigenvalue problem iA*x=lambda*iB*x, using Intlab's function "verifyeig".
%eigs and vects should contain approximate eigenvalues and eigenvectors 
%respectively.

M=length(eigs);
ieigs=intval(NaN(M,1));
if not(all_eigs) %In this case, we only try to validate the first eigenvalue (we might end up validating a couple of them if they are clustered)
    test=1;
    if M>1 && abs((eigs(2)-eigs(1))/eigs(1))<tol_simple_eig
        mult_eig=2;
    else
        mult_eig=1;
    end
    while mult_eig<=M && test
        lambda=sum(eigs(1:mult_eig))/mult_eig;
        [ilambda,~] = verifyeig(iA,lambda,vects(:,1:mult_eig),iB);
        test=max(isnan(ilambda)); 
        if test %The eigenpair(s) could not be validated, trying with a larger cluster
            mult_eig=mult_eig+1;
        end
    end    
    if not(test)        
        ieigs(1:mult_eig)=ilambda;
    else
        error(['Some eigenvalues could not be validated. If this happens, you can either try to increase the coefficient margin '...
               'in homotopy.m, or increase the number of modes K and/or N. You can also try to play with the coefficient inu_r in initialize_tD_L.m'])
    end
else %In this case, we try to validate all the eigenvalues
    m=1;
    while m<=M
        nb_eig_bwd=0;
        nb_eig_fwd=0;
        while m+nb_eig_fwd<M && abs((eigs(m+nb_eig_fwd+1)-eigs(m))/eigs(m))<tol_simple_eig
            nb_eig_fwd=nb_eig_fwd+1;
        end
        test=1;
        while nb_eig_bwd<m && nb_eig_fwd<=M-m && test
            lambda=sum(eigs(m+(-nb_eig_bwd:nb_eig_fwd)))/(1+nb_eig_fwd+nb_eig_bwd);
            [ilambda,~] = verifyeig(iA,lambda,vects(:,m+(-nb_eig_bwd:nb_eig_fwd)),iB);
            test=max(isnan(ilambda));
            if test %The eigenpair(s) could not be validated, trying with a larger cluster
                if m-nb_eig_bwd==1
                    nb_eig_fwd=nb_eig_fwd+1;
                elseif m+nb_eig_fwd==M
                    nb_eig_bwd=nb_eig_bwd+1;
                else
                    if abs((eigs(m-nb_eig_bwd)-eigs(m-nb_eig_bwd-1))) < abs((eigs(m+nb_eig_fwd+1)-eigs(m+nb_eig_fwd)))
                        nb_eig_bwd=nb_eig_bwd+1;
                    else
                        nb_eig_fwd=nb_eig_fwd+1;
                    end
                end
            end
        end
        if not(test)
            ieigs(m+(-nb_eig_bwd:nb_eig_fwd))=ilambda;
            m=m+nb_eig_fwd+1;
        else
            error(['Some eigenvalues could not be validated. If this happens, you can either try to increase the coefficient margin '...
                   'in homotopy.m, or increase the number of modes K and/or N. You can also try to play with the coefficient inu_r in initialize_tD_L.m'])
        end      
    end
end
