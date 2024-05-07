close all;
clear all;

% Define the file paths
t1w_path = 'Linear Regression\subject_1_t1w_aligned_register.nii';
t2w_path = 'Linear Regression\subject_1_t2w_register.nii';
applied_t1w_path = 'Linear Regression\subject_2_t1w_aligned_register.nii';

% Load the NIfTI files
t1w_img = niftiread(t1w_path);
t2w_img = niftiread(t2w_path);
applied_t1w_img = niftiread(applied_t1w_path);

% Extract the 95th z-slice
z_slice = 95;
t1w_slice = t1w_img(:, :, z_slice);
t2w_slice = t2w_img(:, :, z_slice);
applied_slice = applied_t1w_img(:,:,z_slice);

% Remove negative values
t1w_slice(t1w_slice < 0) = 0;
t2w_slice(t2w_slice < 0) = 0;
applied_slice(applied_slice<0) = 0;

% Flatten the 2D slices to 1D arrays for regression
t1w_flat = t1w_slice(:);
t2w_flat = t2w_slice(:);
t1w_flat = applied_slice(:);

% Prepare the data for regression
X = t1w_flat;
y = t2w_flat;

% Create a linear regression model
model = fitlm(X, y);

% Use predict function to automatically apply the linear transformation
% We reshape back to the 2D slice format
predicted_t2w_slice = reshape(predict(model, X), size(applied_slice));

% Plotting the original and predicted slices for comparison
figure;
subplot(1, 3, 1);
imshow(t1w_slice, []);
title('Original T1w Slice');

subplot(1, 3, 2);
imshow(t2w_slice, []);
title('Actual T2w Slice');

subplot(1, 3, 3);
imshow(predicted_t2w_slice, []);
title('Predicted T2w Slice');


% Plotting the regression
figure;
scatter(X, y, 'filled');
hold on;
plot(X, predict(model, X), 'r', 'LineWidth', 2);
xlabel('T1w Intensities');
ylabel('T2w Intensities');
title('Regression Analysis');
legend({'Actual Data', 'Regression Line'}, 'Location', 'best');

% Display regression coefficients
coefficients = table2array(model.Coefficients(:, 'Estimate'));
fprintf('Regression coefficients: Slope = %f, Intercept = %f\n', coefficients(1), coefficients(2));
