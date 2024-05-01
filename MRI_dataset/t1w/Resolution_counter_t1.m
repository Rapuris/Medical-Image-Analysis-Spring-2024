% Directory where your NIfI (.nii) files are stored
directory = 'C:\Users\willi\OneDrive - The Webb Schools\Documents\Sophmore Spring\Med Image Analysis\mia_s24_mri_project\t1w';

% Get a list of .nii files in the directory
files = dir(fullfile(directory, '*.nii.gz'));

% Initialize a map to count files by their x and y resolutions
resolutionCount = containers.Map('KeyType', 'char', 'ValueType', 'int32');

% Traverse through all .nii files
for i = 1:length(files)
    % Construct the full file path
    filePath = fullfile(files(i).folder, files(i).name);
    
    % Read the NIfTI file
    info = niftiinfo(filePath);
   
    % Extract x and y dimensions from the file
    xRes = info.ImageSize(1);
    yRes = info.ImageSize(2);
    zRes = info.ImageSize(3);
    
    % Create a unique key for the resolution
    resolutionKey = sprintf('%dx%dx%d', xRes, yRes, zRes);
    
    % Update the count for this resolution
    if isKey(resolutionCount, resolutionKey)
        resolutionCount(resolutionKey) = resolutionCount(resolutionKey) + 1;
    else
        resolutionCount(resolutionKey) = 1;
    end
end

% Display the count of files for each x and y resolution
keys = resolutionCount.keys;
for i = 1:length(keys)
    fprintf('%s has %d file(s)\n', keys{i}, resolutionCount(keys{i}));
end

disp("_____________________________________________________________");

% Initialize a map to count files by their rounded voxel dimensions
voxelDimCount = containers.Map('KeyType', 'char', 'ValueType', 'int32');

% Traverse through all .nii files
for i = 1:length(files)
    % Construct the full file path
    filePath = fullfile(files(i).folder, files(i).name);
    
    % Read the NIfTI file's header to get voxel dimensions
    info = niftiinfo(filePath);
    
    % Extract voxel dimensions and round to the nearest millimeter
    voxelDims = round(info.PixelDimensions * 10) / 10;
    % voxelDims = info.PixelDimensions;

    % Create a unique key for the voxel dimensions
    voxelKey = sprintf('%0.1fmm x %0.1fmm x %0.1fmm', voxelDims(1), voxelDims(2), voxelDims(3));
    % fprintf('%dmm x %dmm x %dmm\n', voxelDims(1), voxelDims(2), voxelDims(3));
    
    % Update the count for these voxel dimensions
    if isKey(voxelDimCount, voxelKey)
        voxelDimCount(voxelKey) = voxelDimCount(voxelKey) + 1;
    else
        voxelDimCount(voxelKey) = 1;
    end
end

% Display the count of files for each set of voxel dimensions
keys = voxelDimCount.keys;
for i = 1:length(keys)
    fprintf('%s %d\n',keys{i}, voxelDimCount(keys{i}));
end

