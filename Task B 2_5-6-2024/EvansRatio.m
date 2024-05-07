%{
    CODE TO CALCULATE EVANS RATIO

%}

clear all;
close all;

%SEGMENT VENTRICLES

w54 = InnerSkullWidthFinder('subject_54_t1w_reg.nii.gz')      %Call with T1w registered NIfTI file
w55 = InnerSkullWidthFinder('subject_55_t1w_reg.nii.gz')
w56 = InnerSkullWidthFinder('subject_56_t1w_reg.nii.gz')
w57 = InnerSkullWidthFinder('subject_57_t1w_reg.nii.gz')
w58 = InnerSkullWidthFinder('subject_58_t1w_reg.nii.gz')