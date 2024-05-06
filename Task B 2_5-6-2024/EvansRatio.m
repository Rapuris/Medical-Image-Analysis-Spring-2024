%{
    CODE TO CALCULATE EVANS RATIO

%}

clear all;
close all;

%SEGMENT VENTRICLES

maxInnerWidth = InnerSkullWidthFinder('subject_54_t1w_reg.nii.gz')