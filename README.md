# Composite-Intuitionistic-Fuzzy-C-Means-Algorithm

Composite-Intuitionistic-Fuzzy-C-Means-Algorithm
MATLAB implementation for image segmentation experiments based on a fuzzy C-means style clustering framework.
Overview
This repository contains MATLAB code, sample images, ground-truth labels, and ablation scripts for image segmentation experiments. The project includes:
a main clustering implementation;
feature extraction utilities;
evaluation code for segmentation metrics;
example test images and label maps;
ablation scripts for comparative analysis.
Repository Structure
```text
.
├── Source\_Code/
│   ├── demo.m
│   ├── CGFFCM.m
│   ├── CGFFCM\_ablation\*.m
│   ├── FeatureExtractor.m
│   ├── Evaluate.m / EvaluateFast.m
│   ├── run\_ablation\_\*.m
│   ├── \*.jpg
│   ├── class\*.mat
│   └── ablation\_results\_\*.csv
├── Test\_Dataset/
│   ├── Images/
│   └── Ground\_Truth/
├── juleifenxi.m
└── plot\_ellipse.m
```
Requirements
Recommended environment:
MATLAB R2020b or later
Image Processing Toolbox
Statistics and Machine Learning Toolbox
Optional:
Fuzzy Logic Toolbox (used by some auxiliary scripts)
Quick Start
1. Clone the repository
```bash
git clone https://github.com/Lianian1/Composite-Intuitionistic-Fuzzy-C-Means-Algorithm.git
cd Composite-Intuitionistic-Fuzzy-C-Means-Algorithm
```
2. Open MATLAB and enter the source directory
```matlab
cd('Source\_Code')
```
3. Run the demo
```matlab
demo
```
Data Format
Input image
format: `.jpg`
Ground truth
format: `.mat`
label-map size should match the image size
Important Note
The current `demo.m` contains absolute local paths. Before public release or use on another machine, please replace them with relative paths.
Example:
```matlab
Img = imread('101027.jpg');
B = load('class101027.mat');
B = B.class101027;
```
Reproducibility
Some scripts use random initialization.
Fixing the random seed may improve reproducibility.
Ablation scripts are provided for several sample images.
