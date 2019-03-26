function [pxy] = Pxy(x,y)
    % With dot products or no?
    numer = -((x-0.25).^2+(y-0.6).^2)/0.08;
    pxy = exp(numer).*chi(x,y);
end