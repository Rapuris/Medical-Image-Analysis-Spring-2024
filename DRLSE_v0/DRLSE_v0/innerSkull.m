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
    % phi=NeumannBoundCond(phi);
    % [phi_x,phi_y]=gradient(phi);
    % s=sqrt(phi_x.^2 + phi_y.^2);
    % smallNumber=1e-10;  
    % Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
    % Ny=phi_y./(s+smallNumber);
    % curvature=div(Nx,Ny);

    regularTerm=adaptedLevelSet(phi, vx, vy);  % MY FXN

    % diracPhi=Dirac(phi,epsilon);
    % areaTerm=diracPhi.*g; % balloon/pressure force
    % edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
    % phi=phi + timestep*(lambda*edgeTerm + alfa*areaTerm);
    phi=phi + timestep*(mu*regularTerm);
end

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

% Initialize force function f
contour = abs(phi) < 1;

% Alternatively, if the contour is not exactly zero but includes a transition around zero:
% contour = bwperim(phi > 0); % This might be necessary if phi smoothly transitions through zero

% Dilate the contour using an 8-connected structuring element
se = strel('square', 3); % 3x3 square structuring element for 8-connectivity
dilatedContour = imdilate(contour, se);

% Initialize force function f to the minimum value for int32
f = zeros(size(phi));  % Set everything to zero
f(dilatedContour) = 1;  % Set only dilated contour areas to one
    
% Define shrinking factor when dot product is negative
shrink_factor = -1;  % You can adjust this value based on desired speed or sensitivity of contraction
log_response = f.*logsig(dotProduct * shrink_factor);
log_response = imgaussfilt(log_response, 1);   
f = log_response;



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




% minDot = min(dotProduct(:));
% maxDot = max(dotProduct(:));
% midPoint = 0.5;  % Midpoint for zero dot product

% % Scale the dot product to [0, 1] with 0.5 as midpoint for zero
% dotProductImage = (dotProduct + maxDot)/(2*maxDot);  % Normalization
% % dotProductImage = dotProductImage * (1 - 2*midPoint) + midPoint;  % Adjust to include a midpoint of 0.5


% % Display the dot product image
% figure; % Open a new figure window
% imshow(dotProductImage, []); % Show the image
% colormap(jet); % Use a colormap to better visualize the values
% colorbar; % Add a colorbar to indicate scaling
% title('Dot Product Image'); % Title the figure
% % Mask where dot product is negative
% mask = dotProduct > 0;

% % Adjust gradients based on the mask
% vx_masked = vx .* mask;
% vy_masked = vy .* mask;

% % Compute divergence using masked gradients
% f = div(dps .* vx_masked - vx_masked, dps .* vy_masked - vy_masked) + 4 * del2(phi);