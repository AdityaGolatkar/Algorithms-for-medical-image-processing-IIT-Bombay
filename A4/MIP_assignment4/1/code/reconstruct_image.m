function [reconstructed_img] = reconstruct_image(img_data,img_mask,c,ref)
imageNoisy = ifft2(img_data);
figure;
imshow(imageNoisy);
title('Noisy Image');
colormap(gca,parula), colorbar;

S = double(img_mask);
%S = img_mask;
%Phantom image
alpha1 = 0.9998;
alpha2 = 0.9999999;
alpha3 = 0.99999999;
gamma2 = 0.0036;
gamma3 = 0.00008;


if(c == 1)
    [reconstructed_img] = get_optimal_x_quadratic(img_data,S,alpha1);
    get_rrmse(ref,reconstructed_img)
    figure;
    imshow(reconstructed_img);
    title('Reconstructed-Phantom-Quadratic');
    colormap(gca,parula), colorbar;
    %imwrite(reconstructed_img,'phantom_quadratic.jpg');
end
if(c == 2)
    [reconstructed_img] = get_optimal_x_huber(img_data,S,alpha2,gamma2);
    get_rrmse(ref,reconstructed_img)
    figure;
    imshow(reconstructed_img);
    title('Reconstructed-Phantom-Huber');
    colormap(gca,parula), colorbar;
    %imwrite(reconstructed_img,'phantom_huber.jpg');
end
if(c == 3)
    [reconstructed_img] = get_optimal_x_log(img_data,S,alpha3,gamma3);
    get_rrmse(ref,reconstructed_img)
    figure;
    imshow(reconstructed_img);
    title('Reconstructed-Phantom-Log');
    colormap(gca,parula), colorbar;
    %imwrite(reconstructed_img,'phantom_log.jpg');
end

        