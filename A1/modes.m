function [] = modes()

load hands2D.mat
%load ellipses2D.mat
%shapes = pointSets;

%=============================================
%Get the three dimensions of the shapes matrix
%=============================================
dim = size(shapes,1);
N = size(shapes,2);
N_sp = size(shapes,3);

%===============================================
%Get the mean and subtract it from the pointsets
%===============================================
[mean,s] = meanshape();
shapes_ms = bsxfun(@minus,shapes,mean);

%===============================================
%Create the X and Y matrix for PCA type analysis
%===============================================
for i=1:N_sp
    temp = shapes(:,:,i);
    Z(:,i) = temp(:);
end

%==================================
%Generate the two covariance matrix
%==================================

cm = Z*(Z');

%================================
%Get the eigen values and vectors
%================================
[v,d] = eig(cm);

%for q=1:size(d)
%    x(q) = d(q,q);
%end
for i=1:size(v,1)/2
    m(1,i) = v(size(v,1),2*(i-1)+1);
    m(2,i) = v(size(v,1),2*(i-1)+2);
end
size(m);
x = m(1,:);
y = m(2,:);
%scatter(x,y,[],[1,0,0]);
plot(x,y,'o-');

