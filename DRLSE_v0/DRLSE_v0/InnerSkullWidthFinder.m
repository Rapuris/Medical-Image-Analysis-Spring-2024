% This Matlab code demonstrates an edge-based active contour model adapted for NIfTI files
% Original code by Chunming Li, adapted for NIfTI and multi-frame processing\


clear all;
close all;

% Load NIfTI file
filename = 'subject_11_t1w_aligned_register.nii'; % specify your NIfTI file name here
outputFilename = strrep(filename, '.nii', 'inner_skull.nii'); % Modified filename for output


info = niftiinfo(filename); % get info from the original file for output

% Read NIfTI file
Img = niftiread(filename); % reads NIfTI image
Img(Img < 100) = 0;
    

% Prepare matrix to store all contours
contours = zeros(size(Img), 'like', Img); % Initialize with the same size and type as Img


for frame = 1:size(Img, 3)
    Img_frame = double(Img(:, :, frame));
    
    % Parameter setting
    timestep=25;  
    mu=0.2/timestep;  
    iter_inner=5;
    iter_outer=20;
    lambda=5; 
    alfa= 1.5;  
    epsilon=1.5; 
    sigma=1.5;     

    G=fspecial('gaussian', 15, sigma);
    Img_smooth = conv2(Img_frame, G, 'same');  
    [Ix, Iy] = gradient(Img_smooth);
    f = Ix.^2 + Iy.^2;
    g = 1./(1+f);  


    % Initialize the level set function as a binary step function
    c0=2;
    % Assuming Img_frame is your image frame of size 192x224
    [x, y] = meshgrid(1:size(Img_frame, 2), 1:size(Img_frame, 1));

    % Define the rectangle parameters
    rectX = [5 size(Img_frame, 2)-5]; % 5 pixels from both sides
    rectY = [5 size(Img_frame, 1)-5]; % 5 pixels from top and bottom
    
    % Initialize the level set function (LSF)
    initialLSF = c0 * ones(size(Img_frame)); % Assumes c0 is defined elsewhere as the outside value

    % Set values inside the rectangle to -c0
    initialLSF(y >= rectY(1) & y <= rectY(2) & x >= rectX(1) & x <= rectX(2)) = -c0;


    % The initialized phi is now set from the modified initialLSF
    phi = initialLSF;
    % Level set evolution for FINDING OUTER SKULL
    for n=1:iter_outer
        phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, 'double-well');
    
    end

    % Assume phi is your initial level set function

    % Binary image of the zero level set
    binaryImg = phi <= 0;

    % Compute distance from the zero level set
    insideDist = bwdist(~binaryImg);  % Distance outside the zero level set
    outsideDist = bwdist(binaryImg);  % Distance inside the zero level set

    % Create signed distance function
    % Assign negative to inside distances, positive to outside distances
    initialLSF = outsideDist - insideDist;

    % Correct the signs according to the initial phi
    initialLSF(phi < 0) = -insideDist(phi < 0);
    initialLSF(phi > 0) = outsideDist(phi > 0);

    phi = initialLSF;


    % Level set evolution for FINDING INNER SKULL
    for n=1:iter_outer
        phi = innerSkull(phi, Img_smooth, lambda, mu, alfa, epsilon, timestep, iter_inner);
        if mod(frame, 30) == 0
            % Plot the updated contour
            imagesc(Img_frame, [0, 1000]); axis off; axis equal; colormap(gray); hold on;
            hContour = contour(phi, [0, 0], 'r', 'LineWidth', 2);
            % Update title
            title(['Final zero level contour, Frame ', num2str(frame), ', Iteration ', num2str(n)]);
            % Refresh figure window to show updates
        drawnow;
        end
    end
    % Store the final contour of phi into the contours matrix
    % Convert phi to binary where phi <= 0 is considered the contour
    contours(:, :, frame) = abs(phi) < 1;
end

niftiwrite(contours, outputFilename, info);