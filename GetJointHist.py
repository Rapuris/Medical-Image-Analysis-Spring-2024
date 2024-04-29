import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

def load_nifti_image(file_path):
    try:
        image = nib.load(file_path)
        data = image.get_fdata()
        # Round data to nearest integer
        data = np.round(data).astype(int)
        print(f"Data shape: {data.shape}")
        print(f"Data type: {data.dtype}")
        print(f"Min value: {np.min(data)}, Max value: {np.max(data)}")
        if np.all(data == 0):
            print(f"Warning: Data in {file_path} contains only zeros.")
        return data
    except Exception as e:
        print(f"Failed to load {file_path}: {e}")
        return None

def plot_joint_histogram(data1, data2, bins=100, title='Joint Histogram'):
    if data1 is None or data2 is None:
        print("Data is missing, cannot plot histogram.")
        return

    # Filter the data to include only values greater than 1
    mask1 = data1 > 1
    mask2 = data2 > 1
    filtered_data1 = data1[mask1 & mask2]  # Use logical AND to keep same indices in both
    filtered_data2 = data2[mask1 & mask2]

    # Calculate the joint histogram
    hist, x_edges, y_edges = np.histogram2d(filtered_data1.ravel(), filtered_data2.ravel(), bins=bins)

    if np.sum(hist) == 0:
        print("Histogram data contains only zeros or no valid data. Check data scaling and bin settings.")
    else:
        print("Histogram generated successfully.")
    
    # Plotting the histogram with a logarithmic color scale
    plt.figure(figsize=(8, 6))
    plt.imshow(np.log1p(hist).T, origin='lower', cmap='gray', 
               extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
               aspect='auto', interpolation='nearest')
    plt.colorbar(label='Log of Counts + 1')
    plt.xlabel('T1w Intensity Values')
    plt.ylabel('T2w Intensity Values')
    plt.title(title)
    plt.show()

# File paths to the NIfTI files
t1w_file_path = r'MRI_dataset\t1W in MNI\subject_1_t1w_aligned2_register.nii'
t2w_file_path = r'MRI_dataset\t2w in MNI\subject_1_t2w_register.nii'

# Loading the image data from NIfTI files
t1w_data = load_nifti_image(t1w_file_path)
t2w_data = load_nifti_image(t2w_file_path)

# Plot the joint histogram
plot_joint_histogram(t1w_data, t2w_data, bins=256)  # You can adjust the number of bins here



# ___________________________________________________________________________________________

# import nibabel as nib
# import numpy as np
# import matplotlib.pyplot as plt

# def load_nifti_image(file_path):
#     try:
#         image = nib.load(file_path)
#         data = image.get_fdata()
#         # Convert data to integer after rounding
#         data = np.round(data).astype(int)
#         return data
#     except Exception as e:
#         print(f"Failed to load {file_path}: {e}")
#         return None

# def plot_histogram(data, bins=100, title='Histogram of Positive Image Intensities'):
#     if data is None:
#         print("No data loaded, cannot plot histogram.")
#         return

#     # Filter out zero and negative values
#     filtered_data = data[data > 10]

#     # Plot the histogram of filtered data
#     plt.figure(figsize=(10, 6))
#     plt.hist(filtered_data.ravel(), bins=bins, color='blue', alpha=0.7)
#     plt.title(title)
#     plt.xlabel('Intensity Values')
#     plt.ylabel('Frequency')
#     plt.grid(True)
#     plt.show()

# # File path to the NIfTI file
# file_path = r'MRI_dataset\t2w in MNI\subject_1_t2w_register.nii'

# # t1w_file_path = r'MRI_dataset\t1W in MNI\subject_1_t1w_aligned2_register.nii'
# # t2w_file_path = r'MRI_dataset\t2w in MNI\subject_1_t2w_register.nii'

# # Loading the image data from NIfTI file
# image_data = load_nifti_image(file_path)

# # Plot the histogram of the image data
# plot_histogram(image_data, bins=256)  # Adjust bins as needed
