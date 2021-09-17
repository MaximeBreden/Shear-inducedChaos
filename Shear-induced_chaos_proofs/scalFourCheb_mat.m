function Mat=scalFourCheb_mat(K,N)

%Computation of the matrix (bilinear form) corresponding to the L^2 scalar 
%product in the Fourier-Cheybyshev representation

k=0:K-1;
Elem=-2*(1./(bsxfun(@plus,k',k).^2-1)+1./(bsxfun(@minus,k',k).^2-1));
Elem(1:2:K,2:2:K)=0;
Elem(2:2:K,1:2:K)=0;
Elem(1,:)=Elem(1,:)/2;
Elem(:,1)=Elem(:,1)/2;

Mat=zeros(K*(2*N-1),K*(2*N-1));
for n=1:2*N-1
    Mat((n-1)*K+(1:K),(n-1)*K+(1:K))=Elem;
end
Mat=sparse(Mat);