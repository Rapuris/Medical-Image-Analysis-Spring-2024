#@title Select and visualize your structures of interest in 3D by using the dropdown menu and clicking "Run Interact". 
from ipywidgets import widgets
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure
import SimpleITK as sitk
import numpy as np
from skimage.measure import label, regionprops
from ipywidgets import interact, IntSlider
import nibabel

pred_data = nib.load('MRI_dataset/Segmented imgs/aparc.DKTatlas aseg.deep.mgz').get_fdata()

labels = [
'Lateral-Ventricle',
'Inf-Lat-Vent',
'Cerebellum-White-Matter',
'Cerebellum-Cortex',
'Thalamus-Proper',
'Caudate',
'Putamen',
'Pallidum',
'3rd-Ventricle',
'4th-Ventricle',
'Brain-Stem',
'Hippocampus',
'Amygdala',
'CSF',
'Accumbens-area',
'VentralDC',
'choroid-plexus',
'caudalanteriorcingulate',
'caudalmiddlefrontal',
'cuneus',
'entorhinal',
'fusiform',
'inferiorparietal',
'inferiortemporal',
'isthmuscingulate',
'lateraloccipital',
'lateralorbitofrontal',
'lingual',
'medialorbitofrontal',
'middletemporal',
'parahippocampal',
'paracentral',
'parsopercularis',
'parsorbitalis',
'parstriangularis',
'pericalcarine',
'postcentral',
'posteriorcingulate',
'precentral',
'precuneus',
'rostralanteriorcingulate',
'rostralmiddlefrontal',
'superiorfrontal',
'superiorparietal',
'superiortemporal',
'supramarginal',
'transversetemporal',
'insula']

labels_lookup = {
'Left-Lateral-Ventricle': 4,
'Left-Inf-Lat-Vent': 5,
'Left-Cerebellum-White-Matter': 7,
'Left-Cerebellum-Cortex': 8,
'Left-Thalamus-Proper': 10,
'Left-Caudate': 11,
'Left-Putamen': 12,
'Left-Pallidum': 13,
'Left-3rd-Ventricle': 14,
'Left-4th-Ventricle': 15,
'Left-Brain-Stem': 16,
'Left-Hippocampus': 17,
'Left-Amygdala': 18,
'Left-CSF': 24,
'Left-Accumbens-area': 26,
'Left-VentralDC': 28,
'Left-choroid-plexus': 31,
'Right-Lateral-Ventricle': 43,
'Right-Inf-Lat-Vent': 44,
'Right-Cerebellum-White-Matter': 46,
'Right-Cerebellum-Cortex': 47,
'Right-Thalamus-Proper': 49,
'Right-Caudate': 50,
'Right-Putamen': 51,
'Right-Pallidum': 52,
'Right-Hippocampus': 53,
'Right-Amygdala': 54,
'Right-Accumbens-area': 58,
'Right-VentralDC': 60,
'Right-choroid-plexus': 63,
'Right-3rd-Ventricle': 14,
'Right-4th-Ventricle': 15,
'Right-Brain-Stem': 16,
'Right-CSF': 24,
'ctx-lh-caudalanteriorcingulate': 1002,
'ctx-lh-caudalmiddlefrontal': 1003,
'ctx-lh-cuneus': 1005,
'ctx-lh-entorhinal': 1006,
'ctx-lh-fusiform': 1007,
'ctx-lh-inferiorparietal': 1008,
'ctx-lh-inferiortemporal': 1009,
'ctx-lh-isthmuscingulate': 1010,
'ctx-lh-lateraloccipital': 1011,
'ctx-lh-lateralorbitofrontal': 1012,
'ctx-lh-lingual': 1013,
'ctx-lh-medialorbitofrontal': 1014,
'ctx-lh-middletemporal': 1015,
'ctx-lh-parahippocampal': 1016,
'ctx-lh-paracentral': 1017,
'ctx-lh-parsopercularis': 1018,
'ctx-lh-parsorbitalis': 1019,
'ctx-lh-parstriangularis': 1020,
'ctx-lh-pericalcarine': 1021,
'ctx-lh-postcentral': 1022,
'ctx-lh-posteriorcingulate': 1023,
'ctx-lh-precentral': 1024,
'ctx-lh-precuneus': 1025,
'ctx-lh-rostralanteriorcingulate': 1026,
'ctx-lh-rostralmiddlefrontal': 1027,
'ctx-lh-superiorfrontal': 1028,
'ctx-lh-superiorparietal': 1029,
'ctx-lh-superiortemporal': 1030,
'ctx-lh-supramarginal': 1031,
'ctx-lh-transversetemporal': 1034,
'ctx-lh-insula': 1035,
'ctx-rh-caudalanteriorcingulate': 2002,
'ctx-rh-caudalmiddlefrontal': 2003,
'ctx-rh-cuneus': 2005,
'ctx-rh-entorhinal': 2006,
'ctx-rh-fusiform': 2007,
'ctx-rh-inferiorparietal': 2008,
'ctx-rh-inferiortemporal': 2009,
'ctx-rh-isthmuscingulate': 2010,
'ctx-rh-lateraloccipital': 2011,
'ctx-rh-lateralorbitofrontal': 2012,
'ctx-rh-lingual': 2013,
'ctx-rh-medialorbitofrontal': 2014,
'ctx-rh-middletemporal': 2015,
'ctx-rh-parahippocampal': 2016,
'ctx-rh-paracentral': 2017,
'ctx-rh-parsopercularis': 2018,
'ctx-rh-parsorbitalis': 2019,
'ctx-rh-parstriangularis': 2020,
'ctx-rh-pericalcarine': 2021,
'ctx-rh-postcentral': 2022,
'ctx-rh-posteriorcingulate': 2023,
'ctx-rh-precentral': 2024,
'ctx-rh-precuneus': 2025,
'ctx-rh-rostralanteriorcingulate': 2026,
'ctx-rh-rostralmiddlefrontal': 2027,
'ctx-rh-superiorfrontal': 2028,
'ctx-rh-superiorparietal': 2029,
'ctx-rh-superiortemporal': 2030,
'ctx-rh-supramarginal': 2031,
'ctx-rh-transversetemporal': 2034,
'ctx-rh-insula': 2035}

def label_lookups(structure, hemi):
  # determine what to plot
  if structure[0].isupper():
    if hemi == "left":
      label = labels_lookup["Left-" + structure]
    elif hemi == "right":
      label = labels_lookup["Right-" + structure]
    else:
      label = [labels_lookup["Left-" + structure], labels_lookup["Right-" + structure]]
  else:
    if hemi == "left":
      label = labels_lookup["ctx-lh-" + structure]
    elif hemi == "right":
      label = labels_lookup["ctx-rh-" + structure]
    else:
      label = [labels_lookup["ctx-lh-" + structure], labels_lookup["ctx-rh-" + structure]]
  return label

@widgets.interact_manual(
    hemisphere=['left', 'right', 'both'], structure=labels)
def plot_3d_plotly_shape(structure, hemisphere, show_mesh=True, crop=True, grid=True):
  print(structure, hemisphere)
  import plotly.graph_objects as go
  label = label_lookups(structure, hemisphere)
  test_cond = np.in1d(pred_data, label).reshape(pred_data.shape)
  roi = np.where(test_cond, 1, 0)
  vert_p, faces_p, normals_p, values_p = measure.marching_cubes(roi, 0, spacing=(1, 1, 1))

  fig = go.Figure(data=[go.Mesh3d(
        x=vert_p[:, 0],
        y=vert_p[:, 1],
        z=vert_p[:, 2],
        # i, j and k give the vertices of triangles
        # here we represent the 4 triangles of the tetrahedron surface
        i=faces_p[:, 0],
        j=faces_p[:, 1],
        k=faces_p[:, 2],
        name='y',
        showscale=True
    )]
  )

  if show_mesh:
    #plot surface triangulation
    tri_vertices = vert_p[faces_p]
    Xe = []
    Ye = []
    Ze = []
    
    for T in tri_vertices:
      Xe += [T[k%3][0] for k in range(4)] + [ None]
      Ye += [T[k%3][1] for k in range(4)] + [ None]
      Ze += [T[k%3][2] for k in range(4)] + [ None]
       
    fig.add_trace(go.Scatter3d(x=Xe,
                     y=Ye,
                     z=Ze,
                     mode='lines',
                     name='',
                     line=dict(color= 'rgb(40,40,40)', width=0.5)))
  if crop:
    scale_min = np.min(vert_p, axis=0)
    scale_max = np.max(vert_p, axis=0)
  else:
    scale_min = [0, 0, 0]
    scale_max = pred_data.shape
  fig.update_layout(
    scene = dict(aspectratio=dict(x=1, y=1, z=1),
        xaxis = dict(range=[scale_min[0], scale_max[0]], visible=grid),
        yaxis = dict(range=[scale_min[1], scale_max[1]], visible=grid),
        zaxis = dict(range=[scale_min[2], scale_max[2]], visible=grid),) 
    )

  fig.show()

label = label_lookups("Lateral-Ventricle", "both")
test_cond = np.in1d(pred_data, label).reshape(pred_data.shape)
roi = np.where(test_cond, 1, 0)
vert_p, faces_p, normals_p, values_p = measure.marching_cubes(roi, 0, spacing=(1, 1, 1))

#_____________________________________________________________________________________________________________________________________
# Load the image
image_path = '/Users/sampathrapuri/Desktop/JHU/CS_Coursework/Medical Image Analysis/Project2/t1w/subject_10_t1w.nii.gz'  # Adjust path to your image file
image = sitk.ReadImage(image_path)

# Convert to a numpy array for manipulation
array = sitk.GetArrayFromImage(image)

#_____________________________________________________________________________________________________________________________________

# Threshold the image to isolate the skull (adjust threshold based on your image contrast)
# This is just one simple way; more sophisticated methods might be needed based on image quality.
threshold_value = array.max() * 0.495  # Example threshold at 50% of max intensity
skull_mask = array > threshold_value

# You might need morphological operations to clean up the mask
from scipy.ndimage import binary_closing, generate_binary_structure

struct = generate_binary_structure(3, 2)  # 3D structuring element
cleaned_mask = binary_closing(skull_mask, structure=struct)


# Label the components
labeled_mask = label(cleaned_mask)

# Extract properties of labeled regions
regions = regionprops(labeled_mask)

# Find the largest region based on area
largest_region = max(regions, key=lambda r: r.area)

# Create a mask for the largest region
cleaned_mask = labeled_mask == largest_region.label

#_____________________________________________________________________________________________________________________________________

image_array = array
mask_array = cleaned_mask
def show_slice(slice_index):
    plt.figure(figsize=(10, 10))
    plt.imshow(image_array[slice_index], cmap='gray', interpolation='none')
    plt.imshow(np.ma.masked_where(mask_array[slice_index] == False, mask_array[slice_index]), cmap='autumn', alpha=0.5)
    plt.title(f'Slice {slice_index}')
    plt.axis('off')
    plt.show()

interact(show_slice, slice_index=IntSlider(min=0, max=image_array.shape[0] - 1, step=1, value=image_array.shape[0] // 2));
#_____________________________________________________________________________________________________________________________________

file = nibabel.load("/Users/sampathrapuri/Desktop/JHU/CS_Coursework/Medical Image Analysis/Project2/t1w/subject_10_t1w.nii.gz")

x,y,z = np.where(test_cond == 1)

#____________________________________GET WIDTHS _________________________________________________________________________________________________

# # Get the indices where the mask is 1
# y_indices, x_indices, z_indices = np.where(test_cond == 1)
# # Initialize a dictionary to store points by z-index
# points_by_z = {}

# # Group points by their z-index
# for x, y, z in zip(x_indices, y_indices, z_indices):
#     if z not in points_by_z:
#         points_by_z[z] = []
#     points_by_z[z].append((x, y))

# # Initialize variables to track the maximum distance and corresponding points across all slices
# global_max_distance = 0
# global_point1 = None
# global_point2 = None
# global_z = None

# # Calculate distances within each z-slice
# for z, points in points_by_z.items():
#     if len(points) > 1:  # Only calculate if there are at least two points
#         points_array = np.array(points)
#         distances = distance_matrix(points_array, points_array)
#         max_distance = np.max(distances)
#         if max_distance > global_max_distance:
#             max_distance_indices = np.unravel_index(np.argmax(distances), distances.shape)
#             global_max_distance = max_distance
#             global_point1 = points_array[max_distance_indices[0]]
#             global_point2 = points_array[max_distance_indices[1]]
#             global_z = z

# # Output the results
# if global_point1 is not None and global_point2 is not None:
#     print(f"Maximum distance: {global_max_distance} between points {global_point1} and {global_point2} on slice {global_z}")
# else:
#     print("No valid points found across slices with more than one point.")


#_____________________________________________________________________________________________________________________________________

# Get the indices where the mask is 1
y_indices, x_indices, z_indices = np.where(cleaned_mask == True)
# Initialize a dictionary to store points by z-index
points_by_z = {}

# Group points by their z-index
for x, y, z in zip(x_indices, y_indices, z_indices):
    if z not in points_by_z:
        points_by_z[z] = []
    points_by_z[z].append((x, y))

# Initialize variables to track the maximum distance and corresponding points across all slices
global_max_distance = 0
global_point1 = None
global_point2 = None
global_z = None

# Calculate distances within each z-slice
for z, points in points_by_z.items():
    if len(points) > 1:  # Only calculate if there are at least two points
        points_array = np.array(points)
        distances = distance_matrix(points_array, points_array)
        max_distance = np.max(distances)
        if max_distance > global_max_distance:
            max_distance_indices = np.unravel_index(np.argmax(distances), distances.shape)
            global_max_distance = max_distance
            global_point1 = points_array[max_distance_indices[0]]
            global_point2 = points_array[max_distance_indices[1]]
            global_z = z

# Output the results
if global_point1 is not None and global_point2 is not None:
    print(f"Maximum distance: {global_max_distance} between points {global_point1} and {global_point2} on slice {global_z}")
else:
    print("No valid points found across slices with more than one point.")

#_____________________________________________________________________________________________________________________________________
brain_mask = nibabel.load("/Users/sampathrapuri/Desktop/fast/test.nii.gz")
brain_mask.shape
max_distance
pred_data.shape