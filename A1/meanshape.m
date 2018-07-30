function [mean_shape,shapes_1] = meanshape()
load hands2D.mat
%load ellipses2D.mat
%shapes = pointSets;
%=============================================
%Get the three dimensions of the shapes matrix
%=============================================
dim = size(shapes,1);
N = size(shapes,2);
N_sp = size(shapes,3);
shapes_1 = shapes;

for j=1:1
%===================
%Get the mean shape
%===================
mean_shape = mean(shapes_1,3);

%=================================================
%Align each shape with respect to the current mean
%=================================================
for i=1:N_sp
    aligned_shape = align(shapes_1(:,:,i),mean_shape);
    norm_as = sqrt(aligned_shape(:)'*aligned_shape(:));
    aligned_shape = aligned_shape/norm_as;
    shapes_1(:,:,i) = aligned_shape;
end
end
mean_shape = mean(shapes_1,3);
norm_mean = sqrt(mean_shape(:)'*mean_shape(:));
mean_shape = mean_shape/norm_mean;
for i=1:N_sp
    aligned_shape = align(shapes_1(:,:,i),mean_shape);
    norm_as = sqrt(aligned_shape(:)'*aligned_shape(:));
    aligned_shape = aligned_shape/norm_as;
    shapes_1(:,:,i) = aligned_shape;
end
