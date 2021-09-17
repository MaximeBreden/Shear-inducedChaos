function Mat=derCheb_mat(K)

v=2*(0:K-1)';
l=floor(K/2);

Mat=full(spdiags(repmat(v,[1,l]),(1:2:K-1),K,K));
% Mat=spdiags(repmat(v,[1,l]),(1:2:K-1),K,K);