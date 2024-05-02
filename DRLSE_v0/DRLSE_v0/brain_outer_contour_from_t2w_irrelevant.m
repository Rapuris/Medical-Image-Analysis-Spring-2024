% Load the required packages for NIfTI file handling
% addpath(genpath('/path_to_nifti_toolbox/')); % Update this path to your NIfTI toolbox

% Read NIfTI file
filename = '11_t2w_stripped_registered.nii.gz'; % Specify your NIfTI file name here
brainImg = niftiread(filename); % Reads NIfTI image

% Process each slice if the image is 3D
for slice = 1:size(brainImg, 3)
    % Extract one slice
    img_slice = brainImg(:, :, slice);

    % Find non-zero pixels (brain region)
    brain_mask = img_slice > 0;

    % Find the boundary of the brain
    boundary_mask = bwperim(brain_mask);

    % Find the maximum intensity in the brain region
    max_intensity = max(img_slice(brain_mask));

    % Create a new image where only the boundary is set with max intensity
    new_img_slice = zeros(size(img_slice)); % Initialize with zeros
    new_img_slice(boundary_mask) = max_intensity;

    % Store modified slice back to the 3D image
    brainImg(:, :, slice) = new_img_slice;

    if(mod(slice, 20) == 0)
        % Display the result (optional, can be commented out for faster processing)
        figure;
        imshow(new_img_slice, []);
        title(['Brain Outline for Slice ', num2str(slice)]);
        colormap('gray'); % Use gray colormap for better visibility
        colorbar; % Show colorbar to interpret intensity values
    end
end

% Optionally, save the modified NIfTI file
niftiwrite(brainImg, '11_t2w_stripped_outline.nii');