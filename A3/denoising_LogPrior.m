function [] = denoising_LogPrior()

%=========================
%Loading the images
%(1)Noisy image is complex
%(2)Noiseless is double
%(3)sz is size of the img
%=========================
load assignmentImageDenoisingPhantom.mat;
sz = size(imageNoisy);


%=================
%Tuning parameters
%=================
alpha = 0.18;
gamma = 4.1;
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
while err > 0.3

%========================================================
%Compute the derivative of the Potential function
%For the case of quadratic it is merely 
%2*xi + 2*summation_k=1-4_(xi - xik) k are the neighbours
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

%/////////////////////
%Gradient computation
%\\\\\\\\\\\\\\\\\\\\\

g1 = gamma - gamma./(1+ds/gamma);
g2 = gamma - gamma./(1+dr/gamma);
g3 = gamma - gamma./(1+dl/gamma);
g4 = gamma - gamma./(1+du/gamma);
g5 = gamma - gamma./(1+dd/gamma);

grad_prior = g1 + g2 + g3 + g4 + g5;
grad_likelihood = 2*(nl_img-imageNoisy);
grad = (1-alpha)*grad_likelihood+alpha*grad_prior;
nl_prev = nl_img;
nl_img = nl_img - tau*grad;             %Image obtained after subtracting the gradient

%===================
%Compute the new obj
%===================
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

obj1 = gamma*dr - gamma*gamma*log(1+dr/gamma);
obj2 = gamma*dl - gamma*gamma*log(1+dl/gamma);
obj3 = gamma*du - gamma*gamma*log(1+du/gamma);
obj4 = gamma*dd - gamma*gamma*log(1+dd/gamma);
obj5 = gamma*dr - gamma*gamma*log(1+ds/gamma);
obj_total = obj1+obj2+obj3+obj4+obj5;

obj = (1-alpha)*(do(:)'*do(:))+alpha*(sum(obj_total(:)));

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
%err = sqrt(diff1)/(sz(1)*sz(2))
err = sqrt(diff1)/(sqrt(imageNoiseless(:)'*imageNoiseless(:)));
end

end;

err
subplot(1,2,1)
imshow(abs(imageNoisy));
subplot(1,2,2)
imshow(abs(nl_img));