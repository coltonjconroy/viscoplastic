M = 0.0; htol = 1e-3;
for i = 1:nelems
    %...x-degrees of freedom first
    xd   = 2.*(1./ELEM(i).jacobian);
    hj   = dhdx(1,i);
    hL2  = [PHI.L2]*dhdx(:,i);

    vR   = [PHI.b2]*dhdx(:,i);
    vL   = [PHI.b1]*dhdx(:,i);
    if i == 1
        hL  = 0;
        hR  = dhdx(1,i+1);
    elseif i == nelems
        hR  = 0;
        hL  = dhdx(1,i-1);
    else
        hL = dhdx(1,i-1);
        hR = dhdx(1,i+1);
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
        s1 = dhdx(2,i);
        s2 = hR - hj;
        s3 = hj - hL;
        dhdx(2,i) = minmod(s1,s2,s3,xd,M);
        if p > 1
            dhdx(3:p+1,i) = 0.;
        end
    end
end