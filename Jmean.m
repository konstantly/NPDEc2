function result = Jmean(chi, psi, dx)
% A function defining J_n,m,p for the mean of all three Js. It takes as input
% the matrix chi, our estimation of psi and the dx.
    result = zeros(size(psi)); 
    result = (Jpp(chi, psi, dx) + Jpx(chi, psi, dx) + Jxp(chi, psi, dx))/3;
end