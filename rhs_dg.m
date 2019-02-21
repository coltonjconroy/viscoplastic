%-----------------------------
dhdx(:,:) = 0;
% Loop over boundaries (nodes)
%-----------------------------
% left bc
hR = [PHI.b1]*h(:,1,irk);   hL = hR;
fhat_h = hR;
dhdx(:,1)  = dhdx(:,1)  - ELEM(1).jacobian*B.b1*fhat_h;
for i = 2:nnodes-1
    jL  = NODE(i).elems(1);       jR  = NODE(i).elems(2);
    % variables evaluated at boundary int points
    hL = [PHI.b2]*h(:,jL,irk);   hR  = [PHI.b1]*h(:,jR,irk);
    % calc flux 
    fhat_h  = hR;
    dhdx(:,jL)  = dhdx(:,jL)  - ELEM(jL).jacobian*B.b2*fhat_h;
    dhdx(:,jR)  = dhdx(:,jR)  - ELEM(jR).jacobian*B.b1*fhat_h;
end
% right bc
hL     = [PHI.b2]*h(:,nelems,irk);   hR = 0;
fhat_h = hR;
dhdx(:,nelems) = dhdx(:,nelems)-ELEM(nelems).jacobian*B.b2*fhat_h;
%---------------
% Element calcs
%--------------
htol = 1e-12;
for j = 1:nelems
    % evaluate variables at adv int. points
    hj         = [PHI.elem]*h(:,j,irk);
    dhdx(:,j)  = dhdx(:,j) - ELEM(j).jacobian*A*(hj);
    dh_ck  = isnan(dhdx(:,j));
    if max(dh_ck) > 0 || max([PHI.L2]*dhdx(:,j)) > 5000.
        display('dhdx solution is NaN')
        keyboard
    end
end
slope_limit_dhdx
for j = 1:nelems
    hj = [PHI.elem]*h(:,j,irk);
    if p == 0
        if hj < htol && h(:,j-1,irk) >= htol
            ii = ii + 1;
            x_l(ii) = X(j)*L;
            t_length(ii) = t_realc;
        end
    else
        x  = [PSI.xa]*X(ELEM(j).nodes)';
        for i = 1:length(hj)-1
            if real(hj(i)) >= htol && real(hj(i+1)) < htol
                ii = ii + 1;
                x_l(ii) = x(i)*L;
                t_length(ii) = t_realc;
            end
        end
    end
    dhdxj = [PHI.elem]*dhdx(:,j);
    % advective flux evaluated at int. points
    ij  = find(abs(dhdxj) < htol);
    Yhj = Yh(hj,dhdxj);
    ik  = find(Yhj < 0);
    Yhj(ij) = hj(ij); Yhj(ik) = 0;
    Fhj  = fh(Yhj,hj,dhdxj);
    Fhj(ij) = 0;
    if abs(imag(Fhj)) >  0
        keyboard
    end
    % source terms evaluated at int. points
    rhs_h(:,j,irk)  = ELEM(j).jacobian*A*(Fhj);
    h_ck  = isnan(rhs_h(:,j,irk));
    if max(h_ck) > 0 || max([PHI.L2]*rhs_h(:,j,irk)) > 5000.
        display('h solution is NaN')
        keyboard
    end
end
%-----------------------------
% Loop over boundaries (nodes)
%-----------------------------
% left bc
hR    = [PHI.b1]*h(:,1,irk);       hL = hR;
dhdxR = [PHI.b1]*dhdx(:,1);     dhdxL = 0; % ?????
if abs(dhdxL) > htol; 
    YhL = Yh(hL,dhdxL); 
    if YhL < 0; YhL = 0; end;
    fhL = fh(YhL,hL,dhdxL);
else
    YhL = hL;
    fhL = 0;
end
if abs(dhdxR) > htol; 
    YhR = Yh(hR,dhdxR); 
    if YhR < 0; YhR = 0; end;
    fhR = fh(YhR,hR,dhdxR);
else
    YhR = hR;
    fhR = 0; 
end     
if abs(dhdxL) > htol;
    lambda1 = speed(YhL,dhdxL);
else
    lambda1 = 0;
end
if abs(dhdxR) > htol;
    lambda2 = speed(YhR,dhdxR);
else
    lambda2 = 0;
end
lambda  = max(lambda1,lambda2);
fhat_h  = 1/2*(fhL  + fhR  - abs(lambda)*(hR - hL));
rhs_h(:,1,irk)  = rhs_h(:,1,irk)  + ELEM(1).jacobian*B.b1*fhat_h;
for i = 2:nnodes-1
    jL  = NODE(i).elems(1);       jR  = NODE(i).elems(2);
    % variables evaluated at boundary int points
    hL    = [PHI.b2]*h(:,jL,irk);   hR  = [PHI.b1]*h(:,jR,irk);
    dhdxL = [PHI.b2]*dhdx(:,jL);  dhdxR = [PHI.b1]*dhdx(:,jR);
    % flux evaluated at boundary int points
    if abs(dhdxL) > htol; 
        YhL = Yh(hL,dhdxL); 
        if YhL < 0; YhL = 0; end;
        fhL  = fh(YhL,hL,dhdxL); 
    else
        YhL = hL; 
        fhL = 0;
    end
    if abs(dhdxR) > htol; 
        YhR = Yh(hR,dhdxR); 
        if YhR < 0; YhR = 0; end;
        fhR = fh(YhR,hR,dhdxR);
        if abs(imag(fhR)) > 0
            keyboard
        end
    else
        YhR = hR;
        fhR = 0;
    end
    % calc eigenvalues
    if abs(dhdxL) > htol;
        lambda1 = speed(YhL,dhdxL);
    else
        lambda1 = 0;
    end
    if abs(dhdxR) > htol;
        lambda2 = speed(YhR,dhdxR);
    else
        lambda2 = 0;
    end
    lambda  = max(abs(lambda1),abs(lambda2));
    % calc flux (local-Lax Friedrichs)
    fhat_h  = 1/2*(fhL  + fhR  - abs(lambda)*(hR - hL));
    rhs_h(:,jL,irk)  = rhs_h(:,jL,irk)  + ELEM(jL).jacobian*B.b2*fhat_h;
    rhs_h(:,jR,irk)  = rhs_h(:,jR,irk)  + ELEM(jR).jacobian*B.b1*fhat_h;
    if abs(imag(fhat_h)) >  0 
        keyboard
    end
end
% right bc
hL  = [PHI.b2]*h(:,nelems,irk);    hR = 0;
dhdxR = [PHI.b1]*dhdx(:,1);     dhdxL = dhdxR;
if abs(dhdxL) > htol; 
    YhL = Yh(hL,dhdxL); 
    if YhL < 0; YhL = 0; end;
    fhL  = fh(YhL,hL,dhdxL);
else
    YhL = hL; 
    fhL = 0;
end
if abs(dhdxR) > htol;
    YhR = Yh(hR,dhdxR); 
    if YhR < 0; YhR = 0; end;
    fhR = fh(YhR,hR,dhdxR);
else
    YhR = hR; 
    fhR = 0;
end
% calc eigenvalues
if abs(dhdxL) > htol;
    lambda1 = speed(YhL,dhdxL);
else
    lambda1 = 0;
end
if abs(dhdxR) > htol;
    lambda2 = speed(YhR,dhdxR);
else
    lambda2 = 0;
end
lambda  = max(lambda1,lambda2);
fhat_h  = 1/2*(fhL  + fhR  - abs(lambda)*(hR - hL));
rhs_h(:,nelems,irk) = rhs_h(:,nelems,irk)+ELEM(nelems).jacobian*B.b2*fhat_h;




