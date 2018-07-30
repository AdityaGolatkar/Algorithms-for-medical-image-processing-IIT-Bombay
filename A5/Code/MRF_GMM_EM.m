function [] = MRF_GMM_EM(imageData,imageMask)

sz = size(imageData);
mean1 = zeros(3,1);
variance1 = ones(3,1);

K = 3;

%Initial estimates of the Prior.
initial_X = zeros(size(imageData));

masked_pixels = imageData(imageMask == 1);
%size(masked_pixels)
mn1 = mean(masked_pixels(:));
var1 = var(masked_pixels(:));

t1 = (mn1 - 1*sqrt(var1));
t2 = (mn1 + 1*sqrt(var1));

initial_X(imageData < t1) = 1;
initial_X((imageData>=t1)&(imageData<t2)) = 2;
initial_X((imageData>=t2)) = 3;

%Finding the mean and variance for each class
for i =1:3
    ind = find(initial_X == i);
    %size(ind)
    m = mean(imageData(ind));
    mean1(i,1) = m;
    v = var(imageData(ind));
    variance1(i,1) = v;
end
%mean1
%variance1

%THus the initial parameter values have been computed along with the
%initial Prior

mem_ship(:,:,1) = zeros(size(imageData));
mem_ship(:,:,2) = zeros(size(imageData));
mem_ship(:,:,3) = zeros(size(imageData));
beta = 0.8;


%Compute the no of ones in the mask and find their position.
%Store them in mask_ones. ones_len will have the length of mask_ones.
%Now we wiil generate a random number from 1 to ones_len and then update
%that positions label.
mask_ones = find(imageMask == 1);
ones_len = length(mask_ones);
check = ones(ones_len,1);
check_cp = check;
all_cases = [1,2,3];

mean_prev = mean1;
variance_prev = variance1;

curr_labels = initial_X;
temp_labels = zeros(size(curr_labels));


i=0;
e=0;
while(sum(abs(mean1-mean_prev)) > 0.001 || sum(abs(variance1 - variance_prev)) > 0.001 || e==0)
    %i = i+1
    e=e+1;
    %This is the potential function for the cliques.
    %potential= zeros(1,3);
    total_post_prob = 0;
    while (true)
    
        if(sum(check) == 0) 
            break;
        end
        %This is a random number in the range 1 to ones_len.
        x_0 = randi(ones_len);

        if(check(x_0,1) == 1)
             
             check(x_0,1) = 0;
        
             y_c = floor((mask_ones(x_0)-1)/sz(1)) + 1;
             x_c = rem(mask_ones(x_0)-1,sz(1)) + 1;
             potential= zeros(1,3);
             if(imageMask(x_c+1,y_c) == 1)
             
                t1 = repmat(curr_labels(x_c+1,y_c),1,3);
                t2 = beta*abs(sign(all_cases - t1));
                potential = potential + t2;
             
             end
             if(imageMask(x_c,y_c+1) == 1)
             
                t1 = repmat(curr_labels(x_c,y_c+1),1,3);
                t2 = beta*abs(sign(all_cases - t1));
                potential = potential + t2;
             
             end
              if(imageMask(x_c-1,y_c) == 1)
             
                t1 = repmat(curr_labels(x_c-1,y_c),1,3);
                t2 = beta*abs(sign(all_cases - t1));
                potential = potential + t2;
             
              end
              if(imageMask(x_c,y_c-1) == 1)
             
                t1 = repmat(curr_labels(x_c,y_c-1),1,3);
                t2 = beta*abs(sign(all_cases - t1));
                potential = potential + t2;
             
              end
             % potential
             prior_prob = exp(-potential);
             
             g(1,1)=exp(-0.5*(imageData(x_c,y_c)-mean1(1))^2/(variance1(1)))/sqrt(variance1(1)); 
             g(1,2)=exp(-0.5*(imageData(x_c,y_c)-mean1(2))^2/(variance1(2)))/sqrt(variance1(2));
             g(1,3)=exp(-0.5*(imageData(x_c,y_c)-mean1(3))^2/(variance1(3)))/sqrt(variance1(3));
             
             post_prob = prior_prob.*g;
             
             [~,new_label] = max(post_prob);
             temp_labels(x_c,y_c) = new_label;
             mem_ship(x_c,y_c,:) = post_prob/(sum(post_prob));
             total_post_prob = total_post_prob + log(post_prob(new_label)/sum(post_prob));
             
        end
        
    
    end
    curr_labels = temp_labels;
    mean_prev = mean1;
    variance_prev = variance1;
    
    %Updating the Mean and Variance.
    for i=1:3
        temp = imageMask.*mem_ship(:,:,i);
        mem_ship(:,:,i);
        mean1(i) =sum(sum(temp.*imageData))/sum(sum(temp));
        temp_1 = repmat(mean1(i),size(imageData));
        temp_2 = (imageData - temp_1).^2;
        variance1(i) = sum(sum(temp.*temp_2))/sum(sum(temp));
    end
    check = check_cp;
    final_post_prob(1,e) = total_post_prob;
end




x1 = find(curr_labels == 1);
x2 = find(curr_labels == 2);
x3 = find(curr_labels == 3);
size(x1);
size(x2);
size(x3);
imageData1 = zeros(sz);
imageData2 = zeros(sz);
imageData3 = zeros(sz);
imageData1(x1) = 255;
imageData1(x2) = 0;
imageData1(x3) = 0;
imageData2(x1) = 0;
imageData2(x2) = 255;
imageData2(x3) = 0;
imageData3(x1) = 0;
imageData3(x2) = 0;
imageData3(x3) = 255;
%img = [imageData1,imageData2,imageData3];
%figure;
img = cat(3,imageData1,imageData2,imageData3);
figure;
imshow(img);
title(['Label Segmented Image for beta = ',num2str(beta)]);
xlabel('RED = Class1, GREEN = Class2, Blue = Class3');
figure
imshow(imageData);
%figure
%imshow(imageMask);
s = length(final_post_prob);
plot(1:s-4,final_post_prob(5:s));
xlabel('Iteration Number');
ylabel('Total Log Posterior Prob');

[~,mem_labels] = max(mem_ship,[],3);
mem_labels = mem_labels.*imageMask;
x1 = find(mem_labels == 1);
x2 = find(mem_labels == 2);
x3 = find(mem_labels == 3);

imageData1 = zeros(sz);
imageData2 = zeros(sz);
imageData3 = zeros(sz);
imageData1(x1) = 255;
imageData1(x2) = 0;
imageData1(x3) = 0;
imageData2(x1) = 0;
imageData2(x2) = 255;
imageData2(x3) = 0;
imageData3(x1) = 0;
imageData3(x2) = 0;
imageData3(x3) = 255;
%img = [imageData1,imageData2,imageData3];
%figure;
img = cat(3,imageData1,imageData2,imageData3);
figure;
imshow(img);
title(['Membership Segmented Image for beta = ',num2str(beta)]);
xlabel('RED = Class1, GREEN = Class2, Blue = Class3');

disp('Class Means');
mean1




