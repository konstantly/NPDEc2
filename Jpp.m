function result = Jpp(chi, psi, dx)
    num = (chi(3:end,2:end-1)-chi(1:end-2,2:end-1)).*(psi(2:end-1,3:end)-psi(2:end-1,1:end-2))...
        - (chi(2:end-1,3:end)-chi(2:end-1,1:end-2)).*(psi(3:end,2:end-1)-psi(1:end-2,2:end-1));
    result = num/(4*dx.*dx);
end