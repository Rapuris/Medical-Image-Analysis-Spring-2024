import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVR
import joblib  # For saving the model

# Define the file paths
t1w_path = 'Linear Regression\\subject_1_t1w_aligned_register.nii'
t2w_path = 'Linear Regression\\subject_1_t2w_register.nii'
output_path = 'Linear Regression\\predicted_subject_1_t2w_register.nii'
model_path = 'Linear Regression\\svr_model.joblib'  # Path to save the trained model

# Load the NIfTI files
t1w_img = nib.load(t1w_path)
t2w_img = nib.load(t2w_path)

# Extract the data as numpy arrays
t1w_data = t1w_img.get_fdata()
t2w_data = t2w_img.get_fdata()

# Ensure no negative values
t1w_data[t1w_data < 0] = 0
t2w_data[t2w_data < 0] = 0

# Flatten the 3D volumes to 1D arrays
t1w_flat = t1w_data.ravel().reshape(-1, 1)  # Keep as 2D array for scikit-learn
t2w_flat = t2w_data.ravel()

# Create an SVM regression model with a polynomial kernel of degree 3
model = SVR(kernel='rbf', C=100, epsilon=0.1)

# Fit the model to the entire 3D volume
model.fit(t1w_flat, t2w_flat)

# Predict using the model for the entire T1w volume
predicted_t2w_flat = model.predict(t1w_flat)
predicted_t2w_data = predicted_t2w_flat.reshape(t1w_data.shape)

# Save the predicted T2w volume as a new NIfTI file
new_t2w_img = nib.Nifti1Image(predicted_t2w_data, t1w_img.affine, t1w_img.header)
nib.save(new_t2w_img, output_path)

# Save the SVM model
joblib.dump(model, model_path)

# Generate a range of T1w intensities from min to max for plotting
t1w_test = np.linspace(t1w_flat.min(), t1w_flat.max(), 300).reshape(-1, 1)
t2w_pred = model.predict(t1w_test)

# Plotting the regression results
plt.figure(figsize=(10, 6))
plt.scatter(t1w_flat, t2w_flat, color='grey', alpha=0.5, label='Actual Data')
plt.plot(t1w_test, t2w_pred, color='red', linewidth=2, label='SVR Prediction')
plt.xlabel('T1w Pixel Intensities')
plt.ylabel('T2w Pixel Intensities')
plt.title('SVR Prediction of T2w from T1w Intensities')
plt.legend()
plt.show()

print("Model and predicted image saved.")
