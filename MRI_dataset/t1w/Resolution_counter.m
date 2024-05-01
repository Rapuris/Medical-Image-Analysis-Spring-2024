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
    
    % Extract x and y dimensions f_t2rom the file
    xRes = info.ImageSize(1);
    yRes = info.ImageSize(2);
    
    % Create a unique key for the resolution
    resolutionKey = sprintf('%dx%d', xRes, yRes);
    
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
    fprintf('Resolution %s has %d file(s)\n', keys{i}, resolutionCount(keys{i}));
end
