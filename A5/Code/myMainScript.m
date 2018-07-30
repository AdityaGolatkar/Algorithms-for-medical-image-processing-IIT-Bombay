tic;
cd ..
cd Data
load assignmentSegmentBrainGmmEmMrf.mat
cd ..
cd Code
MRF_GMM_EM(imageData,imageMask);
toc;