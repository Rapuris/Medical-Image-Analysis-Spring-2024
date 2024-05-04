function phi = innerSkull(phi_0, g, lambda,mu, alfa, epsilon, timestep, iter)
%  This Matlab code implements an edge-based active contour model as an
%  application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper:
%
%      C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation", 
%        IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
%
%  Input:
%      phi_0: level set function to be updated by level set evolution
%      g: edge indicator function
%      mu: weight of distance regularization term
%      timestep: time step
%      lambda: weight of the weighted length term   (smoothness)
%      alfa:   weight of the weighted area term     (balloon force)
%      epsilon: width of Dirac Delta function
%      iter: number of iterations
%      potentialFunction: choice of potential function in distance regularization term. 
%              As mentioned in the above paper, two choices are provided: potentialFunction='single-well' or
%              potentialFunction='double-well', which correspond to the potential functions p1 (single-well) 
%              and p2 (double-well), respectively.%
%  Output:
%      phi: updated level set function after level set evolution
%
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/

phi=phi_0;
g = imdiffusefilt(g);
[vx, vy]=gradient(g);

% vx = imgaussfilt(vx, 1);
% vy = imgaussfilt(vy, 1);
% % {SHOW GRADIENT MAP %}
% x = 1:size(phi, 2);
% y = 1:size(phi, 1);
% figure
% contour(x,y,x.*g)
% hold on
% quiver(x,y,vx,vy)
% hold off

for k=1:iter
    phi=NeumannBoundCond(phi);
    [phi_x,phi_y]=gradient(phi);
    s=sqrt(phi_x.^2 + phi_y.^2);
    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
    Ny=phi_y./(s+smallNumber);
    curvature=div(Nx,Ny);

    if mod(k, 5) == 0
        % % {SHOW CURVATURE MAP %}
        % x = 1:size(phi, 2);
        % y = 1:size(phi, 1);
        % figure
        % contour(phi_x,phi_y,x.*g)
        % hold on
        % quiver(x,y,phi_x,phi_y)
        % hold off
    end
    regularTerm=adaptedLevelSet(phi, vx, vy);  % MY FXN
    % diracPhi=Dirac(phi,epsilon);
    edgeTerm = findEdgeTerm(phi, curvature);  % MY FXN

    delta_phi = edgeTerm .* (1 - regularTerm);
    phi=phi + timestep*(delta_phi);
    % phi=phi + timestep*(edgeTerm + regularTerm);
end

function f = findEdgeTerm(phi, curvature)
    % Initialize f with zeros, assuming curvature is a matrix
    f = zeros(size(curvature));
    % Create a mask for pixels near the level set
    nearLevelSet = abs(phi) < 1.5;
    % Apply logistic sigmoid to curvature values near the level set
    f(nearLevelSet) = logsig(curvature(nearLevelSet));
    % Apply a Gaussian filter to the entire field
    % Filtering here is done to the whole matrix but mainly affects regions where f was modified
    f = imgaussfilt(f, 1);



function f = adaptedLevelSet(phi, vx, vy)
%{ updated part to set negative gradients to 0%}
center = [96, 112];  % Center of the image
[rows, cols] = size(phi);
[X, Y] = meshgrid(1:cols, 1:rows);

% Vectors from center to each pixel
vecX = X - center(2);
vecY = Y - center(1);

% Check the dot product
dotProduct = vecX .* vx + vecY .* vy;
s = sqrt(vx.^2 + vy.^2);

% Initialize force function f
contour = abs(phi) < 1;

f = zeros(size(phi));  % Set everything to zero
f(contour) = 1;  % Set only dilated contour areas to one
    
% Define shrinking factor when dot product is negative
shrink_factor = 1;  % You can adjust this value based on desired speed or sensitivity of contraction
response = f.*logsig(dotProduct * shrink_factor);
response = imgaussfilt(response, 1);   
f = response;



function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  