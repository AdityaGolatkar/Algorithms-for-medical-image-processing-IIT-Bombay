function [sigma] = get_sigma(img)
% This function returns the ML estimate of the standard devaition of the
% i.i.d Gaussian Noise by using the formula derived in the report.
% We are considering the intensities in the background(i.e. air) where the
% true intensity should be 0.
i = 1:30;
j = 1:30;
R1 = img(i,j);
R1 = R1(:);
R1 = abs(R1);
R1 = R1.*R1;
n1 = size(R1,1);
s1 = sum(R1);

i = 226:256;
j = 1:30;
R2 = img(i,j);
R2 = R2(:);
R2 = abs(R2);
R2 = R2.*R2;
n2 = size(R2,1);
s2 = sum(R2);

i = 1:30;
j = 236:256;
R3 = img(i,j);
R3 = R3(:);
R3 = abs(R3);
R3 = R3.*R3;
n3 = size(R3,1);
s3 = sum(R3);

i = 226:256;
j = 236:256;
R4 = img(i,j);
R4 = R4(:);
R4 = abs(R4);
R4 = R4.*R4;
n4 = size(R4,1);
s4 = sum(R4);

sigma_square = (s1+s2+s3+s4)/(n1+n2+n3+n4);
sigma = sqrt(sigma_square);

