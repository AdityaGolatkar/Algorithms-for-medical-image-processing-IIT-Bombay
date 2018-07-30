function [] = denoising_HuberPrior()

%=========================
%Loading the images
%(1)Noisy image is complex
%(2)Noiseless is double
%(3)sz is size of the img
%========================
load assignmentImageDenoisingPhantom.mat;
sz = size(imageNoisy);
magni = abs(imageNoisy);
x = sum(magni(:))/(sz(1)*sz(2));


%=================
%Tuning parameters
%=================
alpha = 0.28;
gamma = 0.3;
tau = 0.5;


%=======================================
%Assuming a Huber Potential function
%=======================================

nl_img = imageNoisy;    %The initial guess
vec_nl_img = nl_img;    %Vectorize the image.This will be needed for Gradient Descent.


%=========================
%loop for gradient descent
%=========================
obj_prev = Inf;                         %This is the objective function to be minimized.        %
err = Inf;

while err > 1

%========================================================
%Compute the derivative of the Potential function
%========================================================

right_shift = circshift(nl_img,[0,1]);  %[0,1] shifts the image right by 1.
left_shift = circshift(nl_img,[0,-1]);  %[0,-1] shifts the image left by 1.
up_shift = circshift(nl_img,[-1,0]);    %[-1,0] shifts the image up by 1.
down_shift = circshift(nl_img,[1,0]);   %[1,0] shifts the image down by 1.


%/////////////////////////////////////////////////////////////
%These are the matrices of image minus the neighbours
%Like dr means each pixel minus its right neighbour and so on
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%ds = abs(nl_img);
%dr = abs(nl_img - right_shift);
%dl = abs(nl_img - left_shift);
%du = abs(nl_img - up_shift);
%dd = abs(nl_img - down_shift);
%do = abs(nl_img - imageNoisy);

ds = nl_img;
dr = nl_img - right_shift;
dl = nl_img - left_shift;
du = nl_img - up_shift;
dd = nl_img - down_shift;
do = nl_img - imageNoisy;

%///////////////////////////////////////////////////////////
%The differences will be the inputs of the g2 function
%Depending on their magnitude the gradient will be decided
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
gs = gradient_g2(ds,gamma);
gr = gradient_g2(dr,gamma);
gl = gradient_g2(dl,gamma);
gu = gradient_g2(du,gamma);
gd = gradient_g2(dd,gamma);

%/////////////////////
%Gradient computation
%\\\\\\\\\\\\\\\\\\\\\

grad_prior = gs + gr + gl + gd + gu;
grad_likelihood = 2*(nl_img-imageNoisy);
grad = (1-alpha)*grad_likelihood+alpha*grad_prior;
nl_prev = nl_img;
nl_img = nl_img - tau*grad;             %Image obtained after subtracting the gradient

%=====================
%Compute the new obj
%=====================
right_shift = circshift(nl_img,[0,1]);  %[0,1] shifts the image right by 1.
left_shift = circshift(nl_img,[0,-1]);  %[0,-1] shifts the image left by 1.
up_shift = circshift(nl_img,[-1,0]);    %[-1,0] shifts the image up by 1.
down_shift = circshift(nl_img,[1,0]);   %[1,0] shifts the image down by 1.

ds = abs(nl_img);
dr = abs(nl_img - right_shift);
dl = abs(nl_img - left_shift);
du = abs(nl_img - up_shift);
dd = abs(nl_img - down_shift);
do = abs(nl_img - abs(imageNoisy));

obj1 = objective(ds,gamma);
obj2 = objective(dr,gamma);
obj3 = objective(dl,gamma);
obj4 = objective(du,gamma);
obj5 = objective(dd,gamma);

obj = (1-alpha)*(do(:)'*do(:))+alpha*(obj1 + obj2 + obj3 + obj4 + obj5);

if obj <= obj_prev
    tau = tau*1.1;
    obj_prev = obj;                   %When the previous prob is greater than current 
    c=0;
else
    tau = tau*0.5;
    nl_img = nl_prev;
    c=1;
end

%===================
%Computing the error
%===================

if c == 0
diff = abs(nl_img)-imageNoiseless;
diff1 = abs(diff(:))'*abs(diff(:));
%err = sqrt(diff1)/(sz(1)*sz(2));
%err;
err = sqrt(diff1)/(sqrt(imageNoiseless(:)'*imageNoiseless(:)));
end

end;

err
subplot(1,2,1)
imshow(abs(imageNoisy));
subplot(1,2,2)
imshow(abs(nl_img));

%======================================================
%Function to compute the gradient of the Huber function
%======================================================
function [g2_] = gradient_g2(img ,gamma)
img1 = abs(img);
g2_ = gamma*ones(size(img));
index = find(img1 <= gamma);
g2_(index) = img(index);

%==========================================
%Function to compute the objective function
%==========================================
function [obj1] = objective(img,gamma)
obj1 = 0;
index = find(img <= gamma);
obj1 = obj1 + 0.5*img(index)'*img(index);
index = find(img > gamma);
x = gamma*img(index) - 0.5*gamma*gamma;
obj1 = obj1 + sum(x);


