%--------------------------------------
% one-dimensional domain discretization
%--------------------------------------
dx = (xN-x0)/((4)*2^(hi-1));              % dx (element length)
X = x0:dx:xN;                             % Nodes
nelems = length(X)-1; nnodes = length(X); % Number of elements & nodes
% Allocate matricies
ELEM(:).nodes = zeros(1,nelems);          % element node info
ELEM(:).jacobian = zeros(1,nelems);       % element Jacobians
% Jacobian calculations (we perform all element calcs on a "master" element)
for i = 1:nelems
    ELEM(i).nodes    = [i,i+1];           % nodes of each element
    ELEM(i).jacobian = 2/(X(i+1)-X(i));   % Jacobian of each element
end
%----
% Set up node-to-element connectivity for periodic bcs (bcs not periodic in
% this case tho)
%----
NODE(:).elems = zeros(nelems, 2);  % Allocate matrix
NODE(1).elems = [ nelems,1 ];      % First node
for i = 2:nnodes-1
    NODE(i).elems = [ i-1,i ];     % Internal nodes
end
NODE(nnodes).elems = [ nelems,1 ]; % Last node