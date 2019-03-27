function result = Jpp(chi, psi, dx)
% What about chi? Should it be input? 
% should it be just read as existing variable?
    result = zeros(size(psi)); %Do we need this?
    result(2:end-1,2:end-1) = ((chi(3:end,2:end-1)-chi(1:end-2,2:end-1)).*...
        (psi(2:end-1,3:end)-psi(2:end-1,1:end-2)) - (chi(2:end-1,3:end)...
        -chi(2:end-1,1:end-2)).*(psi(3:end,2:end-1)-...
        psi(1:end-2,2:end-1)))/(4*dx.*dx);
end