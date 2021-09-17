function ft=evalCheb(f,t,a,b)
%Normalization: f(t)=f_0 + 2\sum_{k\geq 2} f_k T_k(2*x/(b-a)+(b+a)/(b-a)).

if nargin<3
    a=-1;
    b=1;
end

k=0:length(f)-1;
M=cos(acos(2*t/(b-a)-(b+a)/(b-a))*k);%Very bad way of coding this
f(2:end)=2*f(2:end);
ft=M*f;