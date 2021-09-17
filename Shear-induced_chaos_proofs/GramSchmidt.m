function w=GramSchmidt(v,Mat_scal)

M=size(v,2);
w=0*v;
w_Mscal=0*transpose(v);
for i=1:M
    vi=v(:,i);
    alpha=w_Mscal(1:i-1,:)*vi;
    wi=vi-w(:,1:i-1)*alpha;
    w_Mscali=wi'*Mat_scal;
    norm_wi=real(sqrt(w_Mscali*wi));
    w(:,i)=wi/norm_wi;
    w_Mscal(i,:)=w_Mscali/norm_wi;
end