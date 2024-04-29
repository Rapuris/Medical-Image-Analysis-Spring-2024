import numpy as np
from scipy.ndimage import shift
from scipy.stats import entropy

# # Define the images
# A = np.array([
#     [0, 0, 0, 0, 0],
#     [0, 4, 6, 0, 0],
#     [0, 2, 4, 6, 0],
#     [0, 4, 6, 0, 0],
#     [0, 0, 0, 0, 0]
# ])

# B = np.array([
#     [5, 5, 5, 5, 5],
#     [5, 5, 3, 2, 5],
#     [5, 5, 1, 4, 2],
#     [5, 5, 4, 2, 5],
#     [5, 5, 5, 5, 5]
# ])

# # Function to compute the joint histogram.
# def compute_joint_histogram(img1, img2):
#     max_val = max(img1.max(), img2.max()) + 1
#     hist = np.zeros((max_val, max_val), dtype=int)
#     for i in range(img1.shape[0]):
#         for j in range(img1.shape[1]):
#             hist[img1[i, j], img2[i, j]] += 1
#     return hist

# # Function to compute mutual information.
# def mutual_information(joint_hist):
#     pxy = joint_hist / joint_hist.sum()
#     px = pxy.sum(axis=1)
#     py = pxy.sum(axis=0)
#     px_py = np.outer(px, py)
#     nzs = pxy > 0
#     return (pxy[nzs] * np.log2(pxy[nzs] / px_py[nzs])).sum()

# # Function to calculate the entropy from the joint histogram.
# def joint_entropy(joint_hist):
#     pxy = joint_hist / joint_hist.sum()
#     return -np.sum(pxy * np.log2(pxy, where=pxy>0))

# # Function to compute the least squares cost
# def least_squares_cost(img1, img2):
#     return ((img1 - img2) ** 2).sum()

# # Compute joint entropy and mutual information for shifts.
# shifts = [-2, -1, 0, 1, 2]
# results = []

# for s in shifts:
#     shifted_A = shift(A, (0, s), cval=0)
#     joint_hist = compute_joint_histogram(shifted_A, B)
#     mi = mutual_information(joint_hist)
#     je = joint_entropy(joint_hist)
#     results.append((s, je, mi))

# # Print the results
# for shift_amount, je, mi in results:
#     print(f"Shift: {shift_amount}, Joint Entropy: {je}, Mutual Information: {mi}")

# # Compute least squares cost for the shifts.
# for s in shifts:
#     shifted_A = shift(A, (0, s), cval=0)
#     lsc = least_squares_cost(shifted_A, B)
#     print(f"Shift: {s}, Least Squares Cost: {lsc}")



#____Q5___________________________________________________________________________________

# Two channel image data as numpy arrays
T1 = np.array([[3, 8, 7],
               [5, 2, 7],
               [1, 5, 6]])

T2 = np.array([[9, 4, 5],
               [7, 9, 6],
               [7, 6, 3]])

# Initial class means
mu1 = np.array([7, 5])
mu2 = np.array([3, 8])

# Function to classify pixels and calculate new means
def k_means_classification(T1, T2, mu1, mu2, iterations=3):
    # Create an array to hold class assignments (1 or 2) for each pixel
    classifications = np.zeros(T1.shape, dtype=int)
    new_mu1 = mu1
    new_mu2 = mu2
    
    for iteration in range(iterations):
        # Classify each pixel
        for i in range(T1.shape[0]):
            for j in range(T1.shape[1]):
                # Calculate the distance from each mean
                dist_to_mu1 = np.linalg.norm(np.array([T1[i, j], T2[i, j]]) - new_mu1)
                dist_to_mu2 = np.linalg.norm(np.array([T1[i, j], T2[i, j]]) - new_mu2)
                
                # Assign the pixel to the closest mean
                if dist_to_mu1 < dist_to_mu2:
                    classifications[i, j] = 1
                else:
                    classifications[i, j] = 2
        
        # Calculate new means
        mask_mu1 = (classifications == 1)
        mask_mu2 = (classifications == 2)
        
        if np.any(mask_mu1):
            new_mu1 = np.array([T1[mask_mu1].mean(), T2[mask_mu1].mean()])
        if np.any(mask_mu2):
            new_mu2 = np.array([T1[mask_mu2].mean(), T2[mask_mu2].mean()])
        
        # Record the new means for this iteration
        means_after_iteration = (new_mu1, new_mu2)
    
    return classifications, means_after_iteration

# Perform K-means classification
classifications, final_means = k_means_classification(T1, T2, mu1, mu2)

classifications, final_means

print("Classifications after 3 iterations:")
print(classifications)
print("\nFinal class means:")
print("mu1:", final_means[0])
print("mu2:", final_means[1])

