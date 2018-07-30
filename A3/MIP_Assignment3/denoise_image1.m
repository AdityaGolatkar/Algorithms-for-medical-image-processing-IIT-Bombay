function [] = denoise_image(c)
load assignmentImageDenoisingPhantom.mat;
img = imageNoisy;
denoised_img = img;
%Brain image
%for quadratic iter = 10
%for huber iter = 5
%for log iter = 5
%alpha1 = 0.08;
%alpha2 = 0.2;
%alpha3 = 0.3;
%gamma2 = 0.2;
%gamma3 = 0.4;

%Phantom image
%for quadratic iter = 3(rrmse - 0.2086)
%for huber iter = 12(rrmse - 0.1062)
%for log iter = 5(rrmse - 0.1521)
%Quadratic a = 0.08 iter = 3.
%huber a = 0.1 iter = 14
%log a = 0.15 iter = 6
%
alpha1 = 0.08;
alpha2 = 0.1;
alpha3 = 0.15;
gamma2 = 0.1;
gamma3 = 0.15;

for iter = 1:6
    occured = zeros(size(img));
    ct = 0;
    while 1
        i = randi([1 256],1,1);
        j = randi([1 256],1,1);
        occured(i,j) = occured(i,j)+1;
        if(occured(i,j) == 1)
            ct = ct+1;
            % 1-up,2-down,3-left,4-right
            i1 = i-1;
            j1 = j;
            i2 = i+1;
            j2 = j;
            i3 = i;
            j3 = j-1;
            i4 = i;
            j4 = j+1;

            %Wrap around cases
            if(i == 1)
                i1 = 256;
            end
            if(i == 256)
                i2 = 1;
            end
            if(j == 1)
                j3 = 256;
            end
            if(j == 256)
                j4 = 1;
            end

            y = denoised_img(i,j);
            x_n1 = denoised_img(i1,j1);
            x_n2 = denoised_img(i2,j2);
            x_n3 = denoised_img(i3,j3);
            x_n4 = denoised_img(i4,j4);

            if(c == 1)
                optimal_x = get_optimal_x_quadratic(y,x_n1,x_n2,x_n3,x_n4,alpha1);
                denoised_img(i,j) = optimal_x;
            end
            if(c == 2)
                optimal_x = get_optimal_x_huber(y,x_n1,x_n2,x_n3,x_n4,alpha2,gamma2);
                denoised_img(i,j) = optimal_x;
            end
            if(c == 3)
                optimal_x = get_optimal_x_log(y,x_n1,x_n2,x_n3,x_n4,alpha3,gamma3);
                denoised_img(i,j) = optimal_x;
            end
        end
        
        if(ct == 256*256)
            break;
        end
    end
end
get_rrmse(imageNoiseless,denoised_img)
subplot(1,2,1);
imshow(abs(imageNoisy));
subplot(1,2,2);
imshow(abs(denoised_img));
        