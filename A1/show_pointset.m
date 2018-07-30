function [] = show_pointset()
%load hands2D.mat
load ellipses2D.mat
shapes = pointSets;
%=============================================
%Get the three dimensions of the shapes matrix
%=============================================
dim = size(shapes,1);
N = size(shapes,2);
N_sp = size(shapes,3);

%======================
%Plot all the pointsets
%======================
for i=1:N_sp
    x=shapes(1,:,i);
    y=shapes(2,:,i);
    %scatter(x,y,[],rand(1,3));
    plot(x,y,'o-');
    hold on;
end