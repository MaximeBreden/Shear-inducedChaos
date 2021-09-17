function M=convomat_exp(u,m_output)
%Compute the matrix associated to the convolution product of a vector of 
%exponential Fourier coeffs. That is, u is assumed to be a vector of exp 
%coeffs, and M is such that, for any vector v of exp coeffs, u*v=Mv.

%The assumed normalizations is f(x)=\sum_{k} f_k exp(i*k*x) 

%The second input is optional, and can be used to enforce the size of M. By
%default the size is the one of u.

m=(length(u)+1)/2;

if nargin==2
    u=[zeros(m_output-m,1);u;zeros(m_output-m,1)];
    m=(length(u)+1)/2;
end

M=toeplitz([u(m:2*m-1);zeros(m-1,1)],transpose([flip(u(1:m));zeros(m-1,1)]));

if nargin==2 && m>m_output
    M=M(1+(m-m_output):2*m-1-(m-m_output),1+(m-m_output):2*m-1-(m-m_output));
end
