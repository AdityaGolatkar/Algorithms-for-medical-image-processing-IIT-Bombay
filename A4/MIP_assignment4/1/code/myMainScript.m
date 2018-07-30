%% MyMainScript

tic;
%% Your code here
cd ../data;
load 'assignmentImageDenoisingPhantom.mat';
cd ../code;
% Put c = 1 for quadratic prior, c = 2 for Huber prior and c = 3 for log prior
c = 1;
% The next function also plots the objective function vs. number of
% iterations apart from denoising the image
imageDenoised = denoise_image(imageNoisy,c);
disp('RRMSE between Noiseless and Denoised images are:');
rrmse = get_rrmse(imageNoiseless,imageDenoised)
figure;
imshow(imageNoisy);
title('Noisy Image');
colormap(gca,parula), colorbar;
figure;
imshow(imageDenoised);
title('Denoised Image');
colormap(gca,parula), colorbar;
figure;
imshow(imageNoiseless);
title('Noiseless Image');
colormap(gca,parula), colorbar;
toc;
