%--------------------------------------------------------------------------
% 
% A 1d RKDG Herschel-Bulkley model for viscoplastic dam break
%  
% Written by: Colton J Conroy @ 480 Convent, NY, NY on 2/18/19
%
% Eqns from Balmforth 2007. J. Non-Newtonian Fluid Mech.
%
%--------------------------------------------------------------------------
clear;
%------------
% User input
%------------
video   = 'off';
fig     = 'on'; 
u_fig   = 'off';
S       = 'on';
limiter = 'on';
hi    = 5;         % controls # of elements (even # of elements)
p     = 1;         % DG polynomial approximation
i_start = 0;
ndof  = p + 1;
dt    = 1/(hi^2)*(1/(p+1)); % time step
T     = 60;        % total simulation run time in secs.
if p == 0
    limiter = 'off';
end
L2_pts = p+2;
L   = 0.06;  % length of reservoir
H   = 0.05;  % height of reservoir
rho = 1000;  % density kg/m^3
g   = 9.81;  % gravitational acceleration
Ty  = 30;     % yield strength
nz  = 0.900; % power-law index for z-stress
K   = 120;   % consistency 
Bn  = Ty*L/(rho*g*H^2); 
% left boundary and right boundary coordinates
x0 = 0; xN = 5;
%--------------------
% set-up video writer
%--------------------
if strcmpi(video,'on')
    writerObj = VideoWriter('1d_debris_fail.avi');
    writerObj.Quality = 100;
    open(writerObj);
end 
%-----------
% functions
%-----------
Yh  = @(h,dhdx)   h - Bn./abs(dhdx);
fh  = @(Y,h,dhdx) -(nz.*abs(dhdx).^(1/nz-1).*Y.^(1+1/nz))./((nz+1).*(2*nz+1)).*...
                  ((1+2*nz).*h-nz.*Y).*dhdx; 
speed = @(Y,dhdx) (Y.^(1/nz + 1).*dhdx.*nz.*abs(dhdx).^(1/nz - 1))./(nz + 1);
%-------------
% Create mesh
%------------- 
domain_discretization_1d        
%-------------------------------------------
% Set up DG 1D matricies and basis functions
%-------------------------------------------
[L2,A,B,C,PHI,PSI] = LDG_matrices_1d(p,S,L2_pts);
%----------------------
% RK time stepping data
%----------------------
[alpha,beta,crk] = RKssp(4,8); nrk = length(alpha); % number or stages
%---------------------
% Initialize variables
%---------------------
h      = zeros([p+1,nelems,nrk+1]); % debris thickness
dhdx   = zeros([p+1,nelems]);      % slope of thickness
rhs_h  = zeros([p+1,nelems,nrk]);   % RHS for thickness
%----------------------------------
% Initial Condition (L2 projection)
%----------------------------------
for j = 1:nelems
    x  = [PSI.xa]*X(ELEM(j).nodes)'; % Initial condition at integration points.
    h_i = zeros(length(x),1); 
    for i = 1:length(x)
        if x(i) <= 1
            h_i(i) = 1;
        else
            h_i(i) = 0;
        end
    end
    h(:,j,1)  = L2*h_i;    
end
irk = 0;
if strcmpi(limiter,'on')
    slope_limiter
end
%------------------------
% Main time stepping loop
%------------------------
umax = 1/2;         % for video purposes
NT = floor(T/dt)-1; % # of time steps
t = 0; m = 1; time = 0; plot_n = 0; t_real = 0;
x_l = zeros(NT+1,1); t_length = zeros(NT+1,1); 
ii = 0;
for n = 0:NT
    if time >= i_start
        for irk = 1:nrk % RK stages
            t = time + crk(irk)*dt; % time in each RK stage
            t_realc = t*(L/H)*(K*L/(rho*g*H^2))^(1/nz);
            rhs_dg
            %calc_pressure
            for i = 1:irk
                h(:,:,irk+1)  = h(:,:,irk+1) + alpha(irk,i)*h(:,:,i) ...
                    + beta(irk,i)*dt*rhs_h(:,:,i);
            end
            if strcmpi(limiter,'on')
                slope_limiter
            end
        end
        h(:,:,1)  = h(:,:,nrk+1);   h(:,:,2:end) = 0;  rhs_h(:,:,:) = 0;
        %-------
        % Video
        %-------
        if strcmpi(video,'on')
            video_v2
        end
        time = time + dt;
        t_real = time*(L/H)*(K*L/(rho*g*H^2))^(1/nz);
    else
        time = time + dt;
    end
end
slope_limit_dhdx
%---------------
% video and figs
%---------------
 if strcmpi(video,'on')
    close(writerObj);
 end
 if strcmpi(fig,'on')
     irk = 1;
     plot_1d
 end