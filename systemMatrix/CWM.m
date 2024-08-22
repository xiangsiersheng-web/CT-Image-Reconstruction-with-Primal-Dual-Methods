
%  nx,ny : 2D image size 
%  nb    : Number of Detector Bins 
% na    : Number of Projection Angles

% OUTPUT: Weight Matrix [nb*na, nx*ny]


% function W = CWM(nx,ny,nb,na)
function W = CWM(drho,dtheta,nx,ny)
% Parameter for considering ray-spacing
na = 180/dtheta;
nb = floor(sqrt(nx^2+ny^2))+1;
ray_pix = drho;


% Parameter for masked reconstruction
mask = true(nx,ny); 

% Coordinates for pixel centers of the image
x = [0:nx-1] - (nx-1)/2;
y = (-1)*([0:ny-1] - (ny-1)/2);	

 % Ex) Top-most pixel's coord = x(1),y(1) 
[x,y] = ndgrid(x, y);

% Filter pixels not to be reconstructed
x = x(mask(:));
y = y(mask(:));

% Setting the number of image pixels to be reconstructed.
np = length(x);		

% Construct column vector containing Projection Angles
angle = [0:na-1]'/na * pi;

% Calculate projected pixel centers. (Projected onto x-axis)

r = cos(angle) * x' + sin(angle) * y';
r = r / ray_pix;	% Take account for ray-spacing
r = r + (nb+1)/2;   % Indices for bins start from 1
lbs = floor(r);        % Find left bins of projected piexl centers.

weight_lval = 1 - (r-lbs);		% weight value for left bins
wieght_rval = 1 - weight_lval;  % weight value for right bins

nc = nx * ny;	% the number of columns of W
ri = lbs + [0:na-1]'*nb*ones(1,np);	% row indices for W
ci = find(mask(:))'; 
ci = ci(ones(1,na),:);              % column indices for W

good = lbs(:) >= 1 & lbs(:) < nb;	% within FOV cases
if any(~good), warning 'FOV too small', end

%Construct a sparse weight matrix.
% WL = sparse(ri(good), ci(good), weight_lval(good), nb*na, nc);	% left bin
% WR = sparse(ri(good)+1, ci(good), wieght_rval(good), nb*na, nc); % right bin
WL = sparse(ri(good), ci(good), weight_lval(good), nb*na, nc);	% left bin
WR = sparse(ri(good)+1, ci(good), wieght_rval(good), nb*na, nc); % right bin
W = WL + WR;

