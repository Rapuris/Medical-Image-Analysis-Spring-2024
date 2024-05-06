% This Matlab code demonstrates an edge-based active contour model adapted for NIfTI files
% Original code by Chunming Li, adapted for NIfTI and multi-frame processing\


clear all;
close all;

% Load NIfTI file
filename = 'subject_55_t1w_reg.nii.gz'; % specify your NIfTI file name here

% Read NIfTI file
Img = niftiread(filename); % reads NIfTI image
Img(Img < 100) = 0;
    
Img = Img(:, :, 90); % Extract a single frame for testing
% Process each frame (assuming Img is 3D)


for frame = 1:size(Img, 3)
    Img_frame = double(Img(:, :, frame));
    
    % Parameter setting
    timestep=25;  
    mu=0.2/timestep;  
    iter_inner=5;
    iter_outer=30;
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

    % % Plot the updated contour
    % imagesc(Img_frame, [0, 1000]); axis off; axis equal; colormap(gray); hold on;
    % hContour = contour(phi, [0, 0], 'r', 'LineWidth', 2);
    % % Update title
    % title('Initial level contour');
    % % Refresh figure window to show updates
    % drawnow;

    % Level set evolution for FINDING OUTER SKULL
    for n=1:iter_outer
        phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, 'double-well');
    
    end

    % % {SHOW GRADIENT MAP %}
    % [vx, vy] = gradient(Img_frame);
    % vx = imgaussfilt(vx, 2);
    % vy = imgaussfilt(vy, 2);
    % x = 1:size(phi, 2);
    % y = 1:size(phi, 1);
    % figure
    % contour(x,y,x.*Img_frame)
    % hold on
    % quiver(x,y,vx,vy)
    % colorbar
    % hold off


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

    % % Visualize the signed distance function
    % figure;
    % imagesc(initialLSF);
    % colormap(jet);
    % colorbar;
    % title('Signed Distance Function');
    % axis image;
    % axis off;


    % Level set evolution for FINDING INNER SKULL
    for n=1:iter_outer
        % Visualize the signed distance function
        % if(n == 1 || n == iter_outer || mod(n, 5) == 0)
        %     figure;
        %     imagesc(phi);
        %     colormap(jet);
        %     colorbar;
        %     title('Signed Distance Function');
        %     axis image;
        %     axis off;
        % end

        phi = innerSkull(phi, Img_smooth, timestep, iter_inner, 0.5);
        % Plot the updated contour
        imagesc(Img_frame, [0, 1000]); axis off; axis equal; colormap(gray); hold on;
        hContour = contour(phi, [0, 0], 'r', 'LineWidth', 2);
        % Update title
        title(['Final zero level contour, Frame ', num2str(frame), ', Iteration ', num2str(n)]);
        % Refresh figure window to show updates
        drawnow;
    end

    % imagesc(Img_frame, [0, 1000]); axis off; axis equal; colormap(gray); hold on;
    % contour(phi, [0, 0], 'r', 'LineWidth', 2); % Display the zero level contour in red

    % Extract the y-coordinates where phi is zero for each x-coordinate
    [rows, cols] = size(phi);
    max_width = 0;
    x_with_max_width = 0;

    for x = 1:cols
        % Calculate the sign of phi for each element in the column
        phi_signs = sign(phi(:, x));
        % Find indices where the sign changes between consecutive elements
        sign_changes = find(diff(phi_signs) ~= 0);

        if ~isempty(sign_changes)
            % Calculate the width as the difference between the maximum and minimum indices of sign change
            width = max(sign_changes) - min(sign_changes);
            if width > max_width
                max_width = width;
                x_with_max_width = x;
            end
        end
    end

    x_with_max_width

    % % Plot vertical bar at the x-coordinate of the maximum width
    % hold on;
    % plot([x_with_max_width, x_with_max_width], [1, rows], 'b-', 'LineWidth', 2);

    % title(['Final zero level contour', ', Max Width = ', num2str(max_width)]);
    % drawnow; % Refresh figure window
end

% figure;
% imagesc(Img_frame, [0, 1000]); axis off; axis equal; colormap(gray); hold on;
% contour(phi, [0, 0], 'r', 'LineWidth', 1);
% title(['Final zero level contour, Frame ', num2str(frame), ', Iteration ', num2str(n)]);
% pause(1); % Pause to view figures (optional, can be removed for faster processing)



    % % Assuming Img_frame is your image frame of size 192x224
    % [rows, cols] = size(Img_frame);
    % [x, y] = meshgrid(1:cols, 1:rows);

    % % Define the rectangle parameters
    % rectX = [5 cols-5]; % 5 pixels from both sides
    % rectY = [5 rows-5]; % 5 pixels from top and bottom

    % % Create binary image for the rectangle
    % binaryImg = zeros(size(Img_frame));
    % binaryImg(rectY(1):rectY(2), rectX(1):rectX(2)) = 1;

    % % Compute distance transform inside (negative inside the shape)
    % insideDist = bwdist(~binaryImg);

    % % Compute distance transform outside (positive outside the shape)
    % outsideDist = bwdist(binaryImg);

    % % Create signed distance function (SDF)
    % initialLSF = outsideDist - insideDist;
