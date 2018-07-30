function [aligned_shape] = align(p1,p2)
%===========================================================
%p1 is aligned with respect p2 and returned in aligned_shape
%===========================================================

load hands2D.mat
%load ellipses2D.mat
%shapes = pointSets;
%=============================================
%Get the three dimensions of the shapes matrix
%=============================================
dim = size(shapes,1);
N = size(shapes,2);
N_sp = size(shapes,3);

%=====================================
%Get the two shapes pointset p1 and p2
%=====================================
%p1 = shapes(:,:,1);
%p2 = shapes(:,:,8);

%================================================================
%Compute the mean of each pointset and subtract from the pointset
%================================================================
mean_p1 = mean(p1,2);
mean_p2 = mean(p2,2);
p1_t = p1 - mean_p1;
p2_t = p2 - mean_p2;

%========================================
%Make the norm of the pointset equal to 1
%========================================
norm_p1t = sqrt((p1_t(:)')*p1_t(:));
p1_ts = p1_t/norm_p1t;
norm_p2t = sqrt((p2_t(:)')*p2_t(:));
p2_ts = p2_t/norm_p1t;

%===================================
%Compute the Cross Covariance Matrix
%===================================
ccm = p1_ts*eye(N)*(p2_ts');
[U,S,V] = svd(ccm);
M = eye(dim);
R = V*(U');
if(det(R)==-1)
    M(dim,dim) = -1
    R = V*M*(U');
end
det(R);

aligned_shape = R*p1_ts;

diff = p2_ts-R*p1_ts; 
sdiff = sqrt(diff.*diff);
(sdiff(:)'*sdiff(:))/N;

%x=p2_ts(1,:);
%y=p2_ts(2,:);
%scatter(x,y,[],[1,0,0]);
%hold on;
%x=p1_ts(1,:);
%y=p1_ts(2,:);
%scatter(x,y,[],[0,1,0]);
%hold on;
%x=aligned_shape(1,:);
%y=aligned_shape(2,:);
%scatter(x,y,[],[0,0,1]);
