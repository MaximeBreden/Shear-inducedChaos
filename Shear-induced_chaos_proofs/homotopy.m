function [ilower_eig_s,iupper_eig_s,eigenvect_s,steps_s]=homotopy(s,steps_s,op,approx_r,eigenval_0,eigenvect_0,ilower_eig_0,margin,tol_simple_eig,show_s,show_enclosure)

%Applies the homotopy method (algorithm 1 of the paper) to the operator
%stored in op. approx_r contains rigorous error bounds for the 1/r and 
%1/r^2 terms appearing in the operator. The starting value is given in s 
%(typically s=0), numerical eigenvalues and eigenvvector for this value
%should be given in eigenval_0 and eigenvect_0, and rigorous lower bounds
%for the eigenvalues should be in ilower_eig_0. margin is a parameter (<1)
%which controls how the breakpoints for the homotopy are selected. 
%tol_simple_eig is a small parameter which controls when, in the validation
%procedure using Intlab's verifyeig, close eigenvalues should be validated
%together as a cluster. show_s and show_enclosure should be equal to 0 or
%1, and control how much information about is displayed during the
%computations.

%Checking that the intput are valid, and that there is "some room" for a
%first homotopy step.
if not(eigenval_0(end-1)<ilower_eig_0(end))
    error('Something seems wrong with the starting point of the homotopy')
else
    ind_eig=length(eigenval_0);
    fprintf('\nind_eig = %d\n',ind_eig)
    eigs_lost=0;
    while ind_eig>2 && abs((eigenval_0(ind_eig-1)-ilower_eig_0(ind_eig))/ilower_eig_0(ind_eig))<1.5*margin
            ind_eig=ind_eig-1;    
            eigs_lost=eigs_lost+1;
    end
    if eigs_lost>1
        fprintf('\n%d eigenvalues losts\n',eigs_lost)
        fprintf('\nind_eig = %d\n',ind_eig)
    end
end

eigenval_s=eigenval_0(1:ind_eig);
eigenvect_s=eigenvect_0(:,1:ind_eig);
ilower_eig_s=ilower_eig_0(1:ind_eig);
correc=10^-4;%parameter used to slightly adapt the computation of the homotopy breakpoints (not crucial)
pert=10^-5;%finite difference step-size used in delta_s_predictor to approximate the derivate of eigenvalues with respect to s
delta_s_min=10^-6;
while s<1 && ind_eig>1 %%% HOMOTOPY
    eigenval_s=eigenval_s(1:ind_eig-1);
    eigenvect_s=op.Reset_for_BC_tens*eigenvect_s(:,1:ind_eig-1);
    nu_min=min(eigenval_s);
    nu_min=min(1.1*nu_min,nu_min-10^2);
    
    %% Computation of the homotopy "breakpoint" s_{k+1}
%     %Old version, without predictor, a bit more reliable but slower
%     threshold_reached=0;
%     delta_s=1/256;
%     [eigenval_s,eigenvect_s,s]=homotopy_breakpoint(s,delta_s,nu_min,(1-margin)*ilower_eig_s(ind_eig),margin,ind_eig-1,eigenval_s,eigenvect_s,op,threshold_reached,show_s);
%     %end of the old version

    %New version, with predictor, faster but slightly less reliable (i.e.
    %the computation might stop because no appropriate breakpoint was
    %selected)
    small_nb_eigs=min(ind_eig-1,max(10,floor(ind_eig/100)));
    [delta_s_pred,correc]=delta_s_predictor(s,(1-margin)*ilower_eig_s(ind_eig),margin,op,small_nb_eigs,pert,delta_s_min,correc,show_s);%cheap approximate computation of delta_s = s_{k+1}-s_{k}     
    delta_s=delta_s_pred;
    if show_s
        fprintf('\nPredicted delta_s : %f\n',delta_s)
    end
    [eigenval_s,eigenvect_s,s]=homotopy_breakpoint_check(s,delta_s,nu_min,(1-margin)*ilower_eig_s(ind_eig),margin,ind_eig-1,eigenval_s,eigenvect_s,op,delta_s_min,show_s);%Checking that the computed delta_s is appropriate
    %end of the new version
    
    steps_s=[steps_s s];
    fprintf('\nBreakpoint for the homotopy at s=%f\n',s)
    
    op.Mat_BC_s_tens=s*op.Mat_BC_tens+(1-s)*op.Mat_BC_0_tens;
    eigenvect_s=op.Mat_BC_s_tens*eigenvect_s;
    eigenvect_s=orthogonalize(eigenvect_s,op.Mat_scal);%Making sure that the numerical eigenvectors are close to orthogonal
    
    is=intval(s);
    op.iMat_BC_s_tens=is*op.iMat_BC_tens+(1-is)*op.iMat_BC_0_tens;
    ieigenvect_s=op.iMat_BC_s_tens*intval(op.Reset_for_BC_tens*eigenvect_s);
    ieigenvect_s=op.Pad*ieigenvect_s;    

    %% Rayleigh-Ritz method, upper bounds  
    is=intval(s);
    iHs=is*op.iH+(1-is)*op.iH0;
    iv=ieigenvect_s;
    iHv=iHs*iv;
    
    switch op.type %For the computation of the error terms coming from the truncation in 1/r
        case 'tD'
            iv_C0=op.coef_C0_tens*abs(iv);
            
            %Error coefficients
            icoef_inv_r2=is*op.icoef_inv_r2*iv;
            %C0 norm
            icoef_inv_r2_fullC0=op.coef_C0_tens*abs(icoef_inv_r2); 
            %Final error term
            coef_err_tDu=op.inv_r2_err*icoef_inv_r2_fullC0;

            iE1=transpose(coef_err_tDu)*iv_C0;
            
        case 'S'          
            iu=iv(1:end-1,:);
            iu_C0=op.coef_C0_tens*abs(iu);
            ilambda_abs=abs(iv(end,:));

            %Error coefficients, order 1
            icoef_inv_r2=((1-is)*op.icoef_S0_inv_r2+is*op.icoef_cA_inv_r2)*iu;
            icoef_der_inv_r2=((1-is)*op.icoef_S0_der_inv_r2+is*op.icoef_cA_der_inv_r2)*iu;
            icoef_der2_inv_r2=((1-is)*op.icoef_S0_der2_inv_r2+is*op.icoef_cA_der2_inv_r2)*iu;
            icoef_inv_r=(is*op.icoef_cA_inv_r)*iu;
            icoef_der_inv_r=(is*op.icoef_cA_der_inv_r)*iu;
            icoef_der2_inv_r=(is*op.icoef_cA_der2_inv_r)*iu;

            %Error coefficients, order 2
            icoef_inv_r2_inv_r2=((1-is)*op.icoef_S0_inv_r2_inv_r2+is*op.icoef_cA_inv_r2_inv_r2)*iu;
            icoef_inv_r2_inv_r=(is*op.icoef_cA_inv_r2_inv_r)*iu;
            icoef_inv_r_inv_r=(is*op.icoef_cA_inv_r_inv_r)*iu;
            icoef_inv_r_der_inv_r2=(is*op.icoef_cA_inv_r_der_inv_r2)*iu;
            icoef_inv_r_der_inv_r=(is*op.icoef_cA_inv_r_der_inv_r)*iu;

            %C0 norm, order 1
            icoef_inv_r2_fullC0=op.coef_C0_tens*abs(icoef_inv_r2);    
            icoef_der_inv_r2_fullC0=op.coef_C0_tens*abs(icoef_der_inv_r2);
            icoef_der2_inv_r2_fullC0=op.coef_C0_tens*abs(icoef_der2_inv_r2);
            icoef_inv_r_fullC0=op.coef_C0_tens*abs(icoef_inv_r);    
            icoef_der_inv_r_fullC0=op.coef_C0_tens*abs(icoef_der_inv_r);
            icoef_der2_inv_r_fullC0=op.coef_C0_tens*abs(icoef_der2_inv_r);

            %C0 norm, order 2
            icoef_inv_r2_inv_r2_fullC0=op.coef_C0_tens*abs(icoef_inv_r2_inv_r2);
            icoef_inv_r2_inv_r_fullC0=op.coef_C0_tens*abs(icoef_inv_r2_inv_r);
            icoef_inv_r_inv_r_fullC0=op.coef_C0_tens*abs(icoef_inv_r_inv_r);
            icoef_inv_r_der_inv_r2_fullC0=op.coef_C0_tens*abs(icoef_inv_r_der_inv_r2);
            icoef_inv_r_der_inv_r_fullC0=op.coef_C0_tens*abs(icoef_inv_r_der_inv_r);

            coef_err_Slambda=is*(approx_r.inv_r2_err*op.icoef_Astar_tu_inv_r2_fullC0+approx_r.inv_r_err*op.icoef_Astar_tu_inv_r_fullC0);
            coef_err_Su=approx_r.inv_r2_err*icoef_inv_r2_fullC0+approx_r.der_inv_r2_err*icoef_der_inv_r2_fullC0+approx_r.der2_inv_r2_err*icoef_der2_inv_r2_fullC0...
                        +approx_r.inv_r_err*icoef_inv_r_fullC0+approx_r.der_inv_r_err*icoef_der_inv_r_fullC0+approx_r.der2_inv_r_err*icoef_der2_inv_r_fullC0...
                        +approx_r.inv_r2_err^2*icoef_inv_r2_inv_r2_fullC0+approx_r.inv_r2_err*approx_r.inv_r_err*icoef_inv_r2_inv_r_fullC0+approx_r.inv_r_err^2*icoef_inv_r_inv_r_fullC0...
                        +approx_r.inv_r2_err*approx_r.der_inv_r_err*icoef_inv_r_der_inv_r2_fullC0+approx_r.inv_r_err*approx_r.der_inv_r_err*icoef_inv_r_der_inv_r_fullC0;


            iE1=op.ixi0*(transpose(coef_err_Su)*iu_C0+coef_err_Slambda*(transpose(ilambda_abs)*iu_C0+transpose(iu_C0)*ilambda_abs));
    end
    
    iA0=iv'*op.iMat_scal*iv;
    iA1=iHv'*op.iMat_scal*iv+cintval(0,iE1.sup);
    
    [upper_vect_S_s,upper_eig_S_s]=eig(iA1.mid,iA0.mid);%Numerical upper bounds
    [upper_eig_S_s,ind_upper_eig_S_s]=sort(real(diag(upper_eig_S_s)));
    upper_vect_S_s=upper_vect_S_s(:,ind_upper_eig_S_s);
    all_eigs=1;
    iupper_eig_s=validate_eigenvalues(upper_eig_S_s,upper_vect_S_s,iA1,iA0,tol_simple_eig,all_eigs);
    iupper_eig_s=sort(sup(real(iupper_eig_s)),'ascend');%Rigorous upper bounds
    
    if iupper_eig_s(ind_eig-1)<ilower_eig_s(ind_eig)
        %% Lehmann Maehly method, lower bounds
        nu=(iupper_eig_s(ind_eig-1)+ilower_eig_s(ind_eig))/2;
        inu=intval(nu);
        
        iHnuv=iHv-inu*iv;    
        switch op.type %For the computation of the error terms coming from the truncation in 1/r
            case 'tD'
                iHnuv_C0=op.coef_C0_tens*abs(iHnuv);
                iE2=transpose(coef_err_tDu)*iHnuv_C0+transpose(iHnuv_C0)*coef_err_tDu+(transpose(coef_err_tDu.*iv_C0))*(coef_err_tDu.*iv_C0);

            case 'S'          
                iHnuu_C0=op.coef_C0_tens*abs(iHnuv(1:end-1,:));
                iHnulambda_abs=abs(iHnuv(end,:));

                iE2=op.ixi0*(transpose(coef_err_Su)*iHnuu_C0+coef_err_Slambda*(transpose(ilambda_abs)*iHnuu_C0+transpose(iu_C0)*iHnulambda_abs)...
                             +transpose(iHnuu_C0)*coef_err_Su+coef_err_Slambda*(transpose(iHnuu_C0)*ilambda_abs+transpose(iHnulambda_abs)*iu_C0)...
                             +(transpose(coef_err_Su)*coef_err_Su+op.ixi0*coef_err_Slambda^2).*(transpose(iu_C0)*iu_C0)...
                             +coef_err_Slambda*(transpose(iu_C0.*coef_err_Su)*ilambda_abs+transpose(ilambda_abs)*(iu_C0.*coef_err_Su))...
                             +coef_err_Slambda^2*transpose(ilambda_abs)*ilambda_abs);
        end
        
        iB1=iA1-inu*iA0;
        iB2=iHnuv'*op.iMat_scal*iHnuv+cintval(0,iE2.sup); 
        [lower_vect_S_s,lower_eig_S_s]=eig(iB2.mid,iB1.mid);%Numerical lower bounds
        [lower_eig_S_s,ind_lower_eig_S_s]=sort(real(diag(lower_eig_S_s)),'descend');
        lower_vect_S_s=lower_vect_S_s(:,ind_lower_eig_S_s);
        
        all_eigs=0;
        if s==1
            all_eigs=1;
        end
        ilower_eig_s=validate_eigenvalues(lower_eig_S_s,lower_vect_S_s,iB2,iB1,tol_simple_eig,all_eigs);        
        ilower_eig_s=real(ilower_eig_s);           
        mask=not(isnan(ilower_eig_s));
        if not(ilower_eig_s(mask)<0)
            error('Unable to guarantuee that the eigenvalues are negative')
        end
        ilower_eig_s=inu+flip(ilower_eig_s);
        ilower_eig_s=inf(ilower_eig_s);
        mask=not(isnan(ilower_eig_s));
        ilower_eig_s(mask)=sort(ilower_eig_s(mask),'ascend');%Rigorous lower bounds
        
        if show_enclosure || s==1
            fprintf('\nEnclosure for the eigenvalues\n')
            [ilower_eig_s eigenval_s iupper_eig_s]
        end
    else
        error('We cannot prove that the last eigenvalue is well separated from the next one. This probably means the rigorous bounds are not tight enough. You may try to increase K and/or N.')
    end
    
    if s<1
        ind_eig=ind_eig-1;
        eigs_lost=1;
        while ind_eig>2 && abs((iupper_eig_s(ind_eig-1)-ilower_eig_s(ind_eig))/ilower_eig_s(ind_eig))<1.5*margin && not(isnan(ilower_eig_s(ind_eig-1)))
            ind_eig=ind_eig-1;    
            eigs_lost=eigs_lost+1;
        end
        if ind_eig>2 && abs((iupper_eig_s(ind_eig-1)-ilower_eig_s(ind_eig))/ilower_eig_s(ind_eig))<1.5*margin
            %If we end up here, it means we need to drop more eigenvalues
            %than those for which we had computed rigorous lower bounds, so
            %we compute rigorous lower bounds for all of them (of course we
            %probably only need a few, so this could be somewhat optimized)
            all_eigs=1;
            ilower_eig_s=validate_eigenvalues(lower_eig_S_s,lower_vect_S_s,iB2,iB1,tol_simple_eig,all_eigs);        
            ilower_eig_s=real(ilower_eig_s);            
            mask=not(isnan(ilower_eig_s));
            if not(ilower_eig_s(mask)<0)
                error('Unable to guarantuee that the eigenvalues are negative')
            end    
            ilower_eig_s=inu+flip(ilower_eig_s);
            ilower_eig_s=inf(ilower_eig_s);
            mask=not(isnan(ilower_eig_s));
            ilower_eig_s(mask)=sort(ilower_eig_s(mask),'ascend');%Rigorous lower bounds
            if show_enclosure
                fprintf('\nEnclosure for the eigenvalues (now with all the lower bounds)\n')
                [ilower_eig_s eigenval_s iupper_eig_s]
            end
        end
        while ind_eig>2 && abs((iupper_eig_s(ind_eig-1)-ilower_eig_s(ind_eig))/ilower_eig_s(ind_eig))<1.5*margin
            ind_eig=ind_eig-1;    
            eigs_lost=eigs_lost+1;
        end
        fprintf('\n')
        if eigs_lost>2
            fprintf('\n%d eigenvalues losts\n',eigs_lost)
        end
        fprintf('\nind_eig = %d\n',ind_eig)
    end   
end

steps_s

if s<1
    error('Could not reach the end of the homotopy')
end
