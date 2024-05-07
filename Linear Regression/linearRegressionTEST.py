import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

# Select the 95th z-slice
z_slice = 95
t1w_slice = t1w_data[:, :, z_slice]
t2w_slice = t2w_data[:, :, z_slice]

t1w_slice[t1w_slice < 0] = 0
t2w_slice[t2w_slice < 0] = 0

# Flatten the 2D slices to 1D arrays for regression
t1w_flat = t1w_slice.ravel()
t2w_flat = t2w_slice.ravel()

# Prepare the data for regression including pixel positions
# Generate grid of pixel positions
x_positions, y_positions = np.meshgrid(np.arange(t1w_slice.shape[1]), np.arange(t1w_slice.shape[0]))
x_flat = x_positions.ravel()
y_flat = y_positions.ravel()

# Combine intensity and positions into a single feature matrix
X = np.column_stack((t1w_flat, x_flat, y_flat))
y = t2w_flat

# Create a linear regression model
model = LinearRegression()
model.fit(X, y)

# Predict using the model and reshape back to the image shape
predicted_t2w_flat = model.predict(X)
predicted_t2w_slice = predicted_t2w_flat.reshape(t1w_slice.shape)


# # Plotting regression results in a 3D plot
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x_flat, y_flat, t2w_flat, c='blue', marker='o', alpha=0.5, label='Actual T2w Intensities')
# ax.scatter(x_flat, y_flat, predicted_t2w_flat, color='red', alpha=0.1, label='Predicted T2w Intensities')

# ax.set_xlabel('X Position')
# ax.set_ylabel('Y Position')
# ax.set_zlabel('T2w Intensities')
# ax.set_title('3D Regression of T1w Intensities with Positions to T2w Intensities')
# ax.legend()

# plt.show()


# Plotting the original and predicted slices for comparison
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].imshow(t1w_slice, cmap='gray')
ax[0].set_title('Original T1w Slice')
ax[0].axis('off')

ax[1].imshow(t2w_slice, cmap='gray')
ax[1].set_title('Actual T2w Slice')
ax[1].axis('off')

ax[2].imshow(predicted_t2w_slice, cmap='gray', vmin=np.min(t2w_slice), vmax=np.max(t2w_slice))
ax[2].set_title('Predicted T2w Slice')
ax[2].axis('off')
plt.show()

# Save the predicted T2w slice as a new NIfTI file
new_t2w_slice_img = nib.Nifti1Image(predicted_t2w_slice, t1w_img.affine, t1w_img.header)
nib.save(new_t2w_slice_img, output_path)

# Print regression coefficients (first two are coefficients for the x and y positions)
print("Regression coefficients:", model.coef_)
print("Intercept:", model.intercept_)
