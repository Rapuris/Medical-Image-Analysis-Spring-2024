import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVR
import joblib  # For saving the model

# Load the saved SVM model
model = joblib.load('Linear Regression\\svr_model.joblib')

# Assuming 'new_t1w_data' is your new T1w image loaded similarly to before
new_t1w_flat = new_t1w_data.ravel().reshape(-1, 1)

# Predict using the loaded model
new_predicted_t2w_flat = model.predict(new_t1w_flat)
new_predicted_t2w_data = new_predicted_t2w_flat.reshape(new_t1w_data.shape)

# Save or process the new predicted T2w image as needed
new_predicted_t2w_img = nib.Nifti1Image(new_predicted_t2w_data, new_t1w_img.affine, new_t1w_img.header)
nib.save(new_predicted_t2w_img, 'path_to_new_predicted_t2w_image.nii')
