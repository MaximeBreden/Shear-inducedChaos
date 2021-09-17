function val=margin_for_homotopy(eigenval,proportion,q,tol_mult)

%Heuristic computation of a good threeshold to use in order to decide which
%value of s_{k+1} to take in the homotopy method (see algorithm 1 in the
%paper). In practice, with the notations of the paper, we take s_{k+1} such
%that \lambda^{(s_{k+1})}_{M-(k+1)} \approx (1-margin) * \lambda^{(s_{k})}_{M-k},
%where margin is the output of this function.

if nargin<2
    if nargin<3
        if nargin<4
            tol_mult=10^-8;
        end
        q=0.1;
    end
    proportion=0.1;
end
    
rel_err=abs(eigenval(2:end)-eigenval(1:end-1))./abs(eigenval(2:end));
nb=min(length(eigenval),max(10,floor(proportion*length(eigenval))));
rel_err=rel_err(end-nb+1:end);
rel_err=rel_err(rel_err>tol_mult);

val=min([quantile(rel_err,q),max(rel_err)/10,1/length(eigenval)]);
