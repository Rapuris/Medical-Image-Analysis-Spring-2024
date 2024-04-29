% Define the directory containing MRI images
folderPath = 'MRI_dataset\testout'; % Update the path as needed
outputFolder = 'MRI_dataset\t1w_axial'; % Define the output folder for axial images
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Create the output folder if it does not exist
end
files = dir(fullfile(folderPath, '*.nii.gz')); % Looking for .nii.gz files

% % Loop over each file in the directory
% for idx = 1:length(files)
%     % Construct full file path
%     filePath = fullfile(files(idx).folder, files(idx).name);
% 
%     % Load the MRI image data
%     mriData = niftiread(filePath);
% 
%     % Transform the image from coronal to axial
%     axialData = permute(mriData, [3 1 2]); % Adjust the permutation as needed
%     axialData = flip(axialData, 1); % Flip vertically (adjust based on your data's orientation)
% 
%     % Define the path for saving the reoriented image in the new folder
%     savePath = fullfile(outputFolder, ['axial_' files(idx).name]);
% 
%     % Save the reoriented image
%     niftiwrite(axialData, savePath);
% 
%     % Display the original and transformed image
%     figure;
%     subplot(1,2,1);
%     imagesc(squeeze(mriData(:, :, round(size(mriData,3)/2)))); % Middle coronal slice
%     title(['Original Coronal View: ' files(idx).name]);
%     axis image;
%     colormap gray;
% 
%     subplot(1,2,2);
%     imagesc(squeeze(axialData(:, :, round(size(axialData,3)/2)))); % Middle axial slice
%     title(['Transformed Axial View: ' files(idx).name]);
%     axis image;
%     colormap gray;
% end


%FOR SAGGITAL

% Loop over each file in the directory
for idx = 1:length(files)
    % Construct full file path
    filePath = fullfile(files(idx).folder, files(idx).name);
    
    % Load the MRI image data
    mriData = niftiread(filePath);
    
    % Transform the image from sagittal to axial
    axialData = permute(mriData, [3 2 1]); % Adjust the permutation as needed
    axialData = flip(axialData, 1); % Flip vertically
    axialData = flip(axialData, 2); % Flip horizontally if needed

    % Define the path for saving the reoriented image in the new folder
    savePath = fullfile(outputFolder, ['axial_' files(idx).name]);

    % Save the reoriented image
    niftiwrite(axialData, savePath);
    
    % Display the original and transformed image
    figure;
    subplot(1,2,1);
    imagesc(squeeze(mriData(:,:,round(size(mriData,3)/2)))); % Middle sagittal slice
    title(['Original Sagittal View: ' files(idx).name]);
    axis image;
    colormap gray;

    subplot(1,2,2);
    imagesc(squeeze(axialData(:,:,round(size(axialData,3)/2)))); % Middle axial slice
    title(['Transformed Axial View: ' files(idx).name]);
    axis image;
    colormap gray;
end

