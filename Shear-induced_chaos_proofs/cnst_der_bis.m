function v=cnst_der_bis(nu1,nu2)

if nu1<=nu2
    error('nu1 should be larger than nu2')
end

alpha=nu2/nu1;
x=(8+4*log(alpha)+sqrt((8+4*log(alpha))^2-64*log(alpha)))/(-16*log(alpha));
x=max(x,1);

v=1/nu1*(1+(2*x+1)*2*x*alpha^(2*x));