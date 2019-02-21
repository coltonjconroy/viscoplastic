function [L2,A,B,C,PHI,PSI,DPHI,L3] = LDG_matrices_1d(p,S,l2_pts)

[xe,we] = gauss_jacobi_quad(p+1,0,0);  % Gauss points & weights for basis integration
[xa,wa] = gauss_jacobi_quad(l2_pts,0,0);   % Gauss points & weights for function integration
[x,w] = gauss_jacobi_quad(p+1,0,0);      % Gauss points & weights for the vertical soln
PSI(1).elem = polyval([-1/2 1/2],xe); PSI(2).elem = polyval([1/2 1/2],xe); % PSI is used to calculate
PSI(1).xa   = polyval([-1/2 1/2],xa); PSI(2).xa = polyval([1/2 1/2],xa);   % IC from exact soln
PSI(1).c = polyval([-1/2 1/2],0);PSI(2).c = polyval([1/2 1/2],0);
% Allocate Matricies
PHI(:).b1 = zeros(p+1,1); PHI(:).b2 = zeros(p+1,1); 
PHI(:).elem = zeros(p+1,length(xe)); PHI(:).L2 = zeros(p+1,length(xa)); PHI(:).L3 = zeros(p+1,length(x));
IL2 = zeros(p+1,length(xa)); IA = zeros(p+1,length(xe)); IC = zeros(p+1,length(xa)); IL3 =  zeros(p+1,length(x));
M = zeros(p+1,p+1); DPHI.b1 = zeros(p+1,1); DPHI.b2 = zeros(p+1,1); 
DPHI(:).L2 = zeros(p+1,length(xa));  

for i = 0:p % for each degree of freedom 0 -> p
    P = jacobi_poly(i,0,0); dPdx = polyder(P); % Legendre Polynomials and Derivatives
    PHI(i+1).b1 = polyval(P,-1); PHI(i+1).b2 = polyval(P,1);  % Poly evaluated at boundaries
    PHI(i+1).c = polyval(P,0); 
    DPHI(i+1).b1 = polyval(dPdx,-1); DPHI(i+1).b2 = polyval(dPdx,1);
    DPHI(i+1).L2 = polyval(dPdx,xa);
    DPHI(i+1).L3 = polyval(dPdx,x);
    DPHI(i+1).elem = polyval(dPdx,xe);
    PHI(i+1).elem = polyval(P,xe); % basis evaluated @ Gauss points for basis integration
    PHI(i+1).L2   = polyval(P,xa); % basis evaluated @ Gauss points for function integration
    PHI(i+1).L3   = polyval(P,x); % basis evaluated @ Gauss points for vertical -> horizontal
    PHI(i+1).dphi = polyval(dPdx,xe);
    IL2(i+1,:) = (wa.*PHI(i+1).L2)'; % L2 integration
    IL3(i+1,:) = (w.*PHI(i+1).L3)'; % Integration for mapping vertical -> horizontal
    IA(i+1,:) = (we.*PHI(i+1).dphi)'; % Dphi integration
    IC(i+1,:) = (wa.*PHI(i+1).L2)'; % Flux function integration
    M(i+1,i+1) = 2/(2*i+1); % Mass matrix integration
end
L2   = M\IL2; % L2 projection matrix
L3   = M\IL3; % vertical -> horizontal
A    = M\IA; % Advection matrix
B.b2 = -M\[PHI.b2]'; B.b1 = M\[PHI.b1]';  % Boundary matrix
BB.b2 = -M\[DPHI.b2]'; BB.b1 = M\[DPHI.b1]';
if strcmpi(S,'on')
    C    = M\IC; % source matrix
else
    C = zeros(p+1,1);
end
    


