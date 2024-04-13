# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:44:01 2024

@author: willi
"""

import os
import json

# Define the path to the 'skin_lesion_dataset' directory
dataset_directory = "C:\\Users\\willi\\OneDrive - The Webb Schools\\Documents\\Sophmore Spring\\Med Image Analysis\\skin_lesion_dataset"

# Function to search for 'melanoma' in a JSON file
def search_for_melanoma_in_json(file_path, word):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            if word in json.dumps(data).lower():  # Convert the JSON data to a string and make it lowercase
                return True
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return False

# Function to scan all subfolders and files within the dataset directory
def scan_files(directory, word):
    count = 0;
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.json'):  # Ensure we are only looking at JSON files
                full_path = os.path.join(root, file)
                if search_for_melanoma_in_json(full_path, word):
                    # print(f'{word} found in: {full_path}')
                    count += 1
    return count
                    

# Start scanning the dataset directory
melanomas = scan_files(dataset_directory, "melanoma")
common = scan_files(dataset_directory, "common")
atypical = scan_files(dataset_directory, "atypical")

fully_sym = scan_files(dataset_directory, "fully symmetric")
ax1_sym = scan_files(dataset_directory, "symmetric in 1 axes")
fully_asym = scan_files(dataset_directory, "fully asymmetric")


print (melanomas, common, atypical)
print (fully_sym, ax1_sym, fully_asym)

