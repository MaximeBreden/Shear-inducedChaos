function w=orthogonalize(v,Mat_scal,tol)

%Making sure that the numerical eigenvectors in v are close to orthogonal.
%In practice the only places where we have to do something seem to be when
%we have double eigenvalues.

if nargin<3
    tol=10^-3;
end

M=size(v,2);
Gram=v'*Mat_scal*v;
norm=sqrt(real(diag(Gram)));
norm_diag=spdiags(1./norm,[0],M,M);
v=v*norm_diag;
Gram=norm_diag*Gram*norm_diag;

err=abs(Gram-diag(diag(Gram)));
err_mask=err>tol;
if not( err_mask-diag(diag(err_mask,1),1)-diag(diag(err_mask,-1),-1) == 0 )
    %This should not happen unless we have eigenvalues of multiplicicty > 2
    w=GramSchmidt(v,Mat_scal);
else
    w=v;
    indices=diag(err_mask,1);
   for i=1:M-1
       if indices(i)
           w(:,i+1)=w(:,i+1)-(w(:,i+1)'*Mat_scal*w(:,i))*w(:,i);
       end
   end
end