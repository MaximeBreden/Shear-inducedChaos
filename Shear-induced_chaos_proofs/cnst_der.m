function v=cnst_der(nu1,nu2)

%Computes a constant v such that, if nu_1 > nu_2 > 1 and u\in \ell^1_{\nu_1},
%|| u' ||_{\ell^1_{\nu_2}} <= v * || u ||_{\ell^1_{\nu_1}}

if nu1<=nu2
    error('nu1 should be larger than nu2')
end   

if nu2<exp(-intval(0.5))*nu1
    v=1/nu1+4*nu2^2/(nu2^2-1)*(nu2/nu1)^2;
else
    v=1/nu1+4*nu2^2/(nu2^2-1)*exp(-intval(1))/(2*log(nu1/nu2));
end

v_bis=cnst_der_bis(nu1,nu2);

if nu2==1
    v=v_bis;
else
    v=min(v,v_bis);
end
