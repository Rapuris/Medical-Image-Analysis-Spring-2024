import nibabel as nib
import numpy as np

def rotate_and_mirror_nifti(file_path, output_path):
    # Load the NIfTI file
    img = nib.load(file_path)
    data = img.get_fdata()

    # Rotate the Z-axis by +90 degrees (counter-clockwise)
    data = np.rot90(data, axes=(0, 1))

    # Rotate the new Y-axis (previously X-axis) by -90 degrees (clockwise)
    data = np.rot90(data, k=3, axes=(0, 2))

    # Mirror the image along the new Z-axis
    data = np.flip(data, axis=2)

    # Create a new NIfTI image
    new_img = nib.Nifti1Image(data, img.affine, img.header)

    # Save the new NIfTI image to the specified output path
    new_img.to_filename(output_path)

# Example usage
file_path = r'MRI_dataset\t1w\subject_6_t1w.nii.gz'  # Replace with your actual file path
output_path = 'MRI_dataset/aligned T1w/subject_6_t1w_aligned.nii.gz'  # Path where the transformed file will be saved
rotate_and_mirror_nifti(file_path, output_path)

output_dir = 'MRI_dataset/aligned T1w'
new_name = 'subject_5_t1w_aligned.nii'
