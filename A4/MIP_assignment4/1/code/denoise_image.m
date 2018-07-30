function [denoised_img] = denoise_image(img,c)
denoised_img = img;
%Phantom image
alpha1 = 0.015;
alpha2 = 0.5971;
alpha3 = 0.4562;
gamma2 = 0.0101;
gamma3 = 0.0314;

%these iteration values have been calculated on the basis of convergence of
%the cost function, i.e. the negative log likelihood
%for quadratic iter = 14
%for huber iter = 14
%for log iter = 6
num_iter = 14;
if(c == 2)
    num_iter = 14;
end
if(c == 3)
    num_iter = 6;
end

J = zeros(num_iter,1);

for iter = 1:num_iter
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
                [optimal_x,J_min] = get_optimal_x_quadratic(y,x_n1,x_n2,x_n3,x_n4,alpha1);
                J(iter,1) = J(iter,1) + J_min;
                denoised_img(i,j) = optimal_x;
            end
            if(c == 2)
                [optimal_x,J_min] = get_optimal_x_huber(y,x_n1,x_n2,x_n3,x_n4,alpha2,gamma2);
                J(iter,1) = J(iter,1) + J_min;
                denoised_img(i,j) = optimal_x;
            end
            if(c == 3)
                [optimal_x,J_min] = get_optimal_x_log(y,x_n1,x_n2,x_n3,x_n4,alpha3,gamma3);
                J(iter,1) = J(iter,1) + J_min;
                denoised_img(i,j) = optimal_x;
            end
        end
        
        if(ct == 256*256)
            break;
        end
    end
end

iter = 1:num_iter;
figure;
plot(iter,J,'-o');
xlabel('Number Of Iterations');
ylabel('Negative Log-Likelihood');

        