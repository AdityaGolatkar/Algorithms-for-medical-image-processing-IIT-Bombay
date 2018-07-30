function [] = show_mean()

load hands2D.mat
%load ellipses2D.mat
%shapes = pointSets;
%=============================================
%Get the three dimensions of the shapes matrix
%=============================================
dim = size(shapes,1);
N = size(shapes,2);
N_sp = size(shapes,3);

[m_shape,s_1] = meanshape();
%======================
%Plot all the pointsets
%======================
for i=1:N_sp
    x=s_1(1,:,i);
    y=s_1(2,:,i);
    scatter(x,y,[],[1,1,0]);
    %plot(x,y,'o-');
    hold on;
end
 
    x=m_shape(1,:);
    y=m_shape(2,:);
    %scatter(x,y,[],[1,0,0]);
    plot(x,y,'o-');