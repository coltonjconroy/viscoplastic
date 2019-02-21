M = 0.0; htol = 1e-3;
for i = 1:nelems
    %...x-degrees of freedom first
    xd   = 2.*(1./ELEM(i).jacobian);
    hj   = h(1,i,irk+1);
    hL2  = [PHI.L2]*h(:,i,irk+1);

    vR   = [PHI.b2]*h(:,i,irk+1);
    vL   = [PHI.b1]*h(:,i,irk+1);
    if i == 1
        hL  = 0;
        hR  = h(1,i+1,irk+1);
    elseif i == nelems
        hR  = 0;
        hL  = h(1,i-1,irk+1);
    else
        hL = h(1,i-1,irk+1);
        hR = h(1,i+1,irk+1);
    end
    s1 = vR - hj;                
    s2 = hR - hj;                
    s3 = hj - hL;                
    HR = minmod(s1,s2,s3,xd,M); 
    HR = hj + HR;                
    s1 = hj - vL;                 
    HL = minmod(s1,s2,s3,xd,M);  
    HL = hj - HL;               
    cR = vR - HR;              
    cL = vL - HL;              
    if abs(cR) > 1.0d-16 || abs(cL) > 1.0d-16
        s1 = h(2,i,irk+1);
        s2 = hR - hj;
        s3 = hj - hL;
        h(2,i,irk+1) = minmod(s1,s2,s3,xd,M);
        if p > 1
            h(3:p+1,i,irk+1) = 0.;
        end
    end
end
        