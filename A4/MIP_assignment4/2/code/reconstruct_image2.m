function [reconstructed_img] = reconstruct_image2(img_data,img_mask,c)
imageNoisy = ifft2(img_data);
figure;
imshow(abs(imageNoisy));
title('Noisy Image');
colormap(gca,hot), colorbar;

S = double(img_mask);
%S = img_mask;
%Phantom image
alpha1 = 0.9998;
alpha2 = 0.9999999;
alpha3 = 0.99999999;
gamma2 = 0.000002;
gamma3 = 0.000002;

if(c == 1)
    [reconstructed_img] = get_optimal_x_quadratic(img_data,S,alpha1);
    figure;
    imshow(abs(reconstructed_img));
    title('Reconstructed-Brain-Quadratic');
    colormap(gca,hot), colorbar;
    %imwrite(reconstructed_img,'phantom_quadratic.jpg');
end
if(c == 2)
    [reconstructed_img] = get_optimal_x_huber(img_data,S,alpha2,gamma2);
    figure;
    imshow(abs(reconstructed_img));
    title('Reconstructed-Brain-Huber');
    colormap(gca,hot), colorbar;
    %imwrite(reconstructed_img,'phantom_huber.jpg');
end
if(c == 3)
    [reconstructed_img] = get_optimal_x_log(img_data,S,alpha3,gamma3);
    figure;
    imshow(abs(reconstructed_img));
    title('Reconstructed-Brain-Log');
    colormap(gca,hot), colorbar;
    %imwrite(reconstructed_img,'phantom_log.jpg');
end

        