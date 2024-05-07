close all;
clear all;

% Step 1: Generate synthetic images
Image1 = rand(100, 100);  % Random image
Image2 = 2 * Image1 + randn(100, 100) * 0.1;  % Strongly correlated with some noise
Image3 = rand(100, 100);  % A new image to be transformed

% Vectorize images
X = Image1(:);
Y = Image2(:);
Image3Vec = Image3(:);

% Step 2: Fit linear regression
model = fitlm(X, Y);

% Display regression coefficients
coefficients = model.Coefficients.Estimate;
fprintf('Regression coefficients: Intercept = %f, Slope = %f\n', coefficients(1), coefficients(2));

% Step 3: Apply regression to a third image
predictedImageVec = predict(model, Image3Vec);

% Reshape the predicted image back to 100x100
predictedImage = reshape(predictedImageVec, 100, 100);

% Display all images for comparison
figure;
subplot(2, 2, 1);
imshow(Image1, []);
title('Image 1 (Input for Model)');

subplot(2, 2, 2);
imshow(Image2, []);
title('Image 2 (Target for Model)');

subplot(2, 2, 3);
imshow(Image3, []);
title('Image 3 (New Input)');

subplot(2, 2, 4);
imshow(predictedImage, []);
title('Predicted Image from Image 3');
