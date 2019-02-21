function u = minmod(a1,a2,a3,M,dx)

check = M*dx^2;
a = [a1; a2; a3];

if abs(a1) <= check
    u = a1;
elseif a1 > 0 && a2 > 0 && a3 > 0
    u = min(a);
elseif a1 < 0 && a2 < 0 && a3 < 0
    u = max(a);
else
    u = 0; 
end
    