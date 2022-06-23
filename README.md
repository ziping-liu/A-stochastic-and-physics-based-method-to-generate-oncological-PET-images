# A-stochastic-and-physics-based-method-to-generate-oncological-PET-images

Open-source code for the proposed stocahstic and physics-based method to generate clinically realistic 2-D oncological PET images.

REFERENCE:

Liu, Z., Laforest, R., Mhlanga, J., Fraum, T.J., Itani, M., Dehdashti, F., Siegel, B.A. and Jha, A.K., 2021, February. Observer study-based evaluation of a stochastic and physics-based method to generate oncological PET images. In Medical Imaging 2021: Image Perception, Observer Performance, and Technology Assessment (Vol. 11599, p. 1159905). International Society for Optics and Photonics.

This folder contains the following main scripts:

main_generate_highres_tumors.m
This file contains the MATLAB code to generate high-resolution tumor images with corresponding low-resolution background image obtained from clinical images

main_generate_DL_train_ground_truth.m
This file contains the MATLAB code to generate ground-truth tumor-fraction area maps based on the above-generated high-resolution tumor images

main_generate_simulated_PET_images.c
This file contains the C code to first generate forward projection of high-resolution tumor images and low-resolution background images, respectively, using simulated PET systems. The forward projections were then added with incorporation of Poisson noise to feed into OSEM reconstruction algorithm to generate eventual simulated PET images.
