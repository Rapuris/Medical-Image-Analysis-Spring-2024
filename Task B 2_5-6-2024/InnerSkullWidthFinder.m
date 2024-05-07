%{

IN: T1w registered NIfTI file
OUT: int for max width of inner skull
    niftii file for mask of inner skull space


%}

function max_width = InnerSkullWidthFinder(filename)
    % Load NIfTI file
    % filename = 'stripped.nii'; % specify your NIfTI file name here
    % filename = 'subject_55_t1w_reg.nii.gz';

    % Read NIfTI file
    Img = niftiread(filename); % reads NIfTI image
    info = niftiinfo(filename); % Read the header information for writing later
    Img(Img < 100) = 0;

    all_max_widths = zeros(size(Img, 3), 1); % Initialize array to store max widths for each frame

    zeroLevelSetImg = false(size(Img)); % Initialize binary 3D image for zero level set

    % Convert logical image to numeric type before writing to NIfTI
    numericZeroLevelSetImg = uint8(zeroLevelSetImg);


    % Define the output filename dynamically based on the input filename    
    [path, name, ~] = fileparts(filename);  % Extract path and name without extension
    outputFilename = fullfile(path, [name '_zero_level_set.nii']);  % Append suffix and re-add extension


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
            phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, 2*iter_inner, 'double-well');
        
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
            phi = innerSkull(phi, Img_smooth, timestep, iter_inner, 0.5);
            % if mod(frame, 30) == 0
            %     % Plot the updated contour
            %     imagesc(Img_frame, [0, 1000]); axis off; axis equal; colormap(gray); hold on;
            %     hContour = contour(phi, [0, 0], 'r', 'LineWidth', 2);
            %     % Update title
            %     title(['Final zero level contour, Frame ', num2str(frame), ', Iteration ', num2str(n)]);
            %     % Refresh figure window to show updates
            % drawnow;
            % end
        end

        % Update binary image for zero level set
        zeroLevelSetImg(:, :, frame) = phi <= 0; % Capture the zero level set for this frame

        % if mod(frame, 30) == 0 % Update the display every 30 frames
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

            all_max_widths(frame) = max_width; % Store the maximum width for this frame
        
            % % Plot vertical bar at the x-coordinate of the maximum width
            % hold on;
            % plot([x_with_max_width, x_with_max_width], [1, rows], 'b-', 'LineWidth', 1);
        
            % title(['Final zero level contour, Frame ', num2str(frame), ', Width = ', num2str(max_width)]);
            % drawnow; % Refresh figure window
        % end

    end

    % Write the numeric image to a NIfTI file
    % niftiwrite(numericZeroLevelSetImg, outputFilename, info);  % Use the original header info for compatibility


    max_width = max(all_max_widths);