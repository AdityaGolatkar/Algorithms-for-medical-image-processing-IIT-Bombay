function [] = denoising_QuadraticPrior()

%=========================
%Loading the images
%(1)Noisy image is complex
%(2)Noiseless is double
%(3)sz is size of the img
%========================
load assignmentImageDenoisingPhantom.mat;
sz = size(imageNoisy);


%=================
%Tuning parameters
%=================
alpha = 0.3;
gamma = 0.7;
tau = 0.5;


%=======================================
%Assuming a Quadratic Potential function
%=======================================

nl_img = imageNoisy;    %The initial guess
vec_nl_img = nl_img;    %Vectorize the image.This will be needed for Gradient Descent.


%=========================
%loop for gradient descent
%=========================
obj_prev = Inf;                         %This is the objective function to be minimized.        %
err = Inf;
c=0;
while err > 1

%========================================================
%Compute the derivative of the Potential function
%For the case of quadratic it is merely 
%2*xi + 2*summation_k=1-4_(xi - xik) k are the neighbours
%========================================================

right_shift = circshift(nl_img,[0,1]);  %[0,1] shifts the image right by 1.
left_shift = circshift(nl_img,[0,-1]);  %[0,-1] shifts the image left by 1.
up_shift = circshift(nl_img,[-1,0]);    %[-1,0] shifts the image up by 1.
down_shift = circshift(nl_img,[1,0]);   %[1,0] shifts the image down by 1.

%/////////////////////
%Gradient computation
%\\\\\\\\\\\\\\\\\\\\\

grad_prior = 10*nl_img - 2*(right_shift+left_shift+up_shift+down_shift);
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
do = abs(nl_img - imageNoisy);

obj = (1-alpha)*(do(:)'*do(:))+alpha*(dr(:)'*dr(:)+dl(:)'*dl(:)+du(:)'*du(:)+dd(:)'*dd(:)+ds(:)'*ds(:));

if obj < obj_prev
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
err = sqrt(diff1)/(sqrt(imageNoiseless(:)'*imageNoiseless(:)));
%err = sqrt(diff1)/(sz(1)*sz(2));
end


end;
err
subplot(1,2,1)
imshow(abs(imageNoisy));
subplot(1,2,2)
imshow((nl_img));