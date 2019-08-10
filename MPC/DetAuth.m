%% DetAuth simulation
% p1 wants to send to p2 a msg
% p3,p4,.. actively corrupted
% pn,pn-1,... omission corrupted

n = 14;
ta = 6;
tw = 2;
tf = 0;
check = true;
3*ta+2*tw+tf
check = 3*ta+2*tw+tf < n

pix = 0;

A = zeros(1,n);
%A(1,:) = logical(randi(2,1,n)-1);

%ran = randi(2,1,n)-1
%logical(ran)
%~logical(ran)

x = pix*ones(1,n);
x(3:3+ta-1) = abs(x(3:3+ta-1)-1);
A = x;
% 5 denotes perp
x(n-tw+1:end) = 5;
A = [A;x]

%
xp = x;
% 7 denotes "n/v"
xp(xp==5) = 7;
A = [A;xp]

%

A3 = A(3,:);
perps = sum(A3(:) == 7);
assoi = sum(A3(:) == 1);
midenika = sum(A3(:) == 0);
if( perps > ta+tw+tf)
    disp('zombie')
elseif(assoi > ta)
    disp('1')
elseif(midenika > ta)
    disp('0')
else
    disp('perp')
end

