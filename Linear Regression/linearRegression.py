import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# Define the file paths
t1w_path = 'Linear Regression\\subject_1_t1w_aligned_register.nii'
t2w_path = 'Linear Regression\\subject_1_t2w_register.nii'
output_path = 'Linear Regression\\predicted_subject_1_t2w_register.nii'

# Load the NIfTI files
t1w_img = nib.load(t1w_path)
t2w_img = nib.load(t2w_path)

# Extract the data as numpy arrays
t1w_data = t1w_img.get_fdata()
t2w_data = t2w_img.get_fdata()

t1w_data[t1w_data < 0] = 0
t2w_data[t2w_data < 0] = 0

z_slice = 95
t1w_slice = t1w_data[:, :, z_slice]
t2w_slice = t2w_data[:, :, z_slice]


# Flatten the 3D images to 1D arrays
t1w_flat = t1w_data.ravel()
t2w_flat = t2w_data.ravel()

# Prepare the data for regression
X = t1w_flat.reshape(-1, 1)
y = t2w_flat

# Create a linear regression model
model = LinearRegression()
model.fit(X, y)


# Manually apply the linear transformation to each pixel in the T1w image
predicted_t2w = np.empty_like(t1w_data)  # Initialize an empty array of the same shape as t1w_data
for i in range(t1w_data.shape[0]):
    for j in range(t1w_data.shape[1]):
        for k in range(t1w_data.shape[2]):
            # Apply linear model manually: y = mx + b
            predicted_t2w[i, j, k] = model.coef_[0] * t1w_data[i, j, k] + model.intercept_

# Reshape the predictions back to the original shape
predicted_t2w = predicted_t2w.reshape(t1w_data.shape)
predicted_t2w[predicted_t2w < 0] = 0

# Save the predicted T2w image as a new NIfTI file
new_t2w_img = nib.Nifti1Image(predicted_t2w, t1w_img.affine, t1w_img.header)
newt2_data = new_t2w_img.get_fdata()
nib.save(new_t2w_img, output_path)

# Plotting the regression
plt.figure(figsize=(10, 6))
plt.scatter(X, y, alpha=0.5, label='Actual Data', s=1)
plt.plot(X.flatten(), predicted_t2w.flatten(), color='red', label='Regression Line')
plt.xlabel('T1w Intensities')
plt.ylabel('T2w Intensities')
plt.title('Linear Regression of T1w to T2w Intensities')
plt.legend()
plt.show()

print("Regression coefficients:", model.coef_)
print("Intercept:", model.intercept_)