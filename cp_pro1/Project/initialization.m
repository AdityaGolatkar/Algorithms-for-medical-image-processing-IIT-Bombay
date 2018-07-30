function [] = grab_cut(inp,gamma,f,gmm_comps)
thresh1 = 15;
sz = size(inp);
thresh = 0.0001;
old_energy = inf;
energy = zeros(1,thresh1);
breakin = 0;
%==============
%INITIALIZATION
%==============

%===================================================================
%Initialize the Trimap i.e. select the background region(not actual
%background but a pseudo version of it.
%===================================================================

imshow(inp);
rec = getrect;
rec = int32(rec);

x1 = max(rec(1),1);
x2 = min(x1+rec(3),sz(2));
y1 = max(rec(2),1);
y2 = min(y1+rec(4),sz(1));


%=============================================================
%Set the alpha for background to zero and unkown region to one
%=============================================================

alpha = zeros(sz(1),sz(2));
alpha(y1:y2,x1:x2) = 1;
alpha = alpha(:);
inp_crop = double(inp(y1:y2,x1:x2,:));

%================================
%Vectorize the initial RGB image
%================================
temp = inp(:);
N = prod(sz(1:2));
inp_vec(:,1) = temp(1:N);
inp_vec(:,2) = temp(N+1:2*N);
inp_vec(:,3) = temp(2*N+1:3*N);

%========================================================================
%Now Extracting the Unknown and the Background Region from the vectorized
%image
%========================================================================

T_F = double(inp_vec(alpha==1,:));
T_B = double(inp_vec(alpha==0,:));


%=====================================================================
%Since we have computed the cropped image. Lets compute the beta value
%=====================================================================

sztest = size(inp_crop);

temp1 = circshift(inp_crop,[-1,0]);
temp2 = circshift(inp_crop,[0,-1]);

temp3 = bsxfun(@minus,temp1,inp_crop);
temp3 = temp3.*temp3;
temp31 = temp3(1:sztest(1)-1,1:sztest(2)-1,:);
temp4 = bsxfun(@minus,temp2,inp_crop);
temp4 = temp4.*temp4;
temp41 = temp4(1:sztest(1)-1,1:sztest(2)-1,:);

temp5 = sum(sum(sum(temp31)))+sum(sum(sum(temp41)));
d1 = sztest(1)-1;
d2 = sztest(2)-1;
d3 = d1*d2*2+d1+d2;
d4 = temp5/d3;
beat1 = 1/(2*d4);
beta = beat1;


%==========================================================================
%Lets initialize the GMM
%As suggested by the paper we take GMM with 5 gaussians for the Foreground
%and the Background
%==========================================================================


back_k = kmeans(T_B,gmm_comps,'Distance','cityblock','Replicates',5);
fore_k = kmeans(T_F,gmm_comps,'Distance','cityblock','Replicates',5);

%==========================================================================
%Using k means we have clustered the foreground and background regions.Each
%cluster will correspond to a gaussian in the Fg and Bg. Who's parameters
%will be sample mean and sample variance as given in SECTION 3.2 of the
%paper
%==========================================================================
pie_b = zeros(1,gmm_comps);
mu_b = zeros(1,3,gmm_comps);
sigma_b = zeros(3,3,gmm_comps);

pie_f = pie_b;
mu_f = mu_b;
sigma_f = sigma_b;

for i=1:gmm_comps
    temp1 = T_B(back_k==i,:);
    mu_b(:,:,i) = double(mean(temp1,1));
    sigma_b(:,:,i) = double(cov(temp1));
    pie_b(i) = size(temp1,1)/size(T_B,1);
    
    temp2 = T_F(fore_k==i,:);
    mu_f(:,:,i) = double(mean(temp2,1));
    sigma_f(:,:,i) = double(cov(temp2));
    pie_f(i) = size(temp2,1)/size(T_F,1);
end


%===========================================================
%Finding all the edges in the Graph along with their weights
%===========================================================
branches = compute_branches(inp_crop,gamma,beta,temp3,temp4);

%===============================================
%Current foreground region is where alpha is one
%===============================================
fore = alpha==1;
back = ~fore;
o = 1;

%=============================
%ITERATIVE ENERGY MINIMIZATION
%=============================
while true
    o
    %======================================================================
    %STEP 1.
    %Now lets assign new K values to each pixel in the image
    %For that we need to compute D for each pixel and each possible Gaussian
    %and choose the minimum
    %======================================================================
   
    
    T_F = inp_vec(fore,:);
    T_B = inp_vec(back,:);
    
    %======================================================================
    %D stores the D values for each pixel for each gaussian component.
    %So we take one gaussian at a time and compute the its D for all pixels
    %The choose the gaussian for each pixel which gives the minimum D.
    %======================================================================
    
    D = zeros(sum(fore),gmm_comps);
    for i=1:gmm_comps
        D(:,i) = -log(mvnpdf(double(T_F),double(mu_f(:,:,i)),double(sigma_f(:,:,i)+f*eye(3)))) -log(pie_f(i));
    end
    [~,fore_k] = min(D,[],2);
    
    D = zeros(sum(back),gmm_comps);
    for i=1:gmm_comps
        D(:,i) = -log(mvnpdf(double(T_B),double(mu_b(:,:,i)),double(sigma_b(:,:,i)+f*eye(3))))-log(pie_b(i))-1.5*log(2*pi);
    end
    [~,back_k] = min(D,[],2);
    
    %We have thus assigned GMM component to each pixel
    
    
    %======================================================================
    %STEP 2.
    %Now for each gaussian in the foreground and background compute the
    %ideal parameter values using sample mean and sample covariance as
    %given in SECTION 3.2 of the paper
    %======================================================================
    
    for i=1:gmm_comps
        temp1 = double(T_B(back_k==i,:));
        mu_b(:,:,i) = double(mean(temp1,1));
        sigma_b(:,:,i) = double(cov(temp1));
        pie_b(i) = size(temp1,1)/size(T_B,1);
        
        temp2 = double(T_F(fore_k==i,:));
        mu_f(:,:,i) = double(mean(temp2,1));
        sigma_f(:,:,i) = double(cov(temp2));
        pie_f(i) = size(temp2,1)/size(T_F,1);
    end
    
    
    %=================
    %STEP 3.
    %Segment the image 
    %=================
    
    sz1 = size(inp_crop);
    no = prod(sz1(1:2));
    
    D_pix = zeros(2,no);
    for y=1:sz1(1)
        for x=1:sz1(2)
            
            pixel(1,1) = inp_crop(y,x,1);
            pixel(1,2) = inp_crop(y,x,2);
            pixel(1,3) = inp_crop(y,x,3);
            pos = (x-1)*sz1(1)+y;
            
            
            for i=1:gmm_comps
                if( isnan(((sigma_f(:,:,i)+sigma_f(:,:,i)')/2)+f*eye(3)) )
                    breakin=1;
                    break;
                end
                F(i) = -log(mvnpdf(pixel,mu_f(:,:,i),(sigma_f(:,:,i)+sigma_f(:,:,i)')/2+f*eye(3)))-log(pie_f(i));
                D_pix(1,pos) = min(F);
                if(isnan((sigma_b(:,:,i)+sigma_b(:,:,i)')/2+f*eye(3)))
                    breakin=1;
                    break;
                end
                F1(i) = -log(mvnpdf(pixel,mu_b(:,:,i),(sigma_b(:,:,i)+sigma_b(:,:,i)')/2+f*eye(3)))-log(pie_b(i));
                D_pix(2,pos) = min(F1);
            end
            if breakin == 1
                break;
            end
        end
            if breakin == 1
                break;
            end
    end
    
    %================================================
    %From here on we will be using the matlab wrapper
    %function for BOkoy Kolmogrov alogorithm
    %================================================
    
    %======================================================
    %The BK_Create will create a Random Field i.e our graph
    %======================================================
    mrf = BK_Create(no);
    
    %================================================================
    %Next we will setup the cost of the individual nodes in the graph
    %================================================================
    BK_SetUnary(mrf,D_pix);
    
    %=================================================================
    %Next we will be assigning the graph edges the potential function
    %values 
    %==================================================================
    BK_SetPairwise(mrf,branches);
    
    %==================================================================
    %Finally the optimal labelling is computed using graph cut and the
    %energy for the labelling is returned
    %==================================================================
    e = BK_Minimize(mrf);
    
    %==========================================================
    %The optimal labelling is obtained using the function below
    %==========================================================
    labels = BK_GetLabeling(mrf);
    
    %===================================
    %To clear up memory delete the graph
    %===================================
    BK_Delete(mrf);
    
    alpha1 = alpha;
    alpha1(alpha1==1) = labels;
    fore(alpha1==2) = 0;
    
    energy(o) = e;
    delta = (old_energy-e)/old_energy;
    
    if breakin == 1
        break;
    end
    
    if delta < thresh
        break;
    end
    
    if o == thresh1
        break;
    end
    
    back = ~fore;
    old_energy = e;
    o=o+1;
    res = inp;
    res_vec = inp_vec;
    res_vec(back,:) = 255;

    for i=1:size(res,2)
        res(:,i,:) = res_vec((i-1)*sz(1)+1:i*sz(1),:);
    end
    imshow(res);
    drawnow;
end

figure
plot(1:o,energy(1:o));
xlabel('Iterations');
ylabel('Energy');
title('Energy Minimzation');
imwrite(res,'man.jpg');

%==========================================
%This function computes the edge potentials
%We are considering 4 neighbourhood system
%==========================================

function branches = compute_branches(im_sub, gamma,betta,temp3,temp4)
temp3 = sum(temp3,3);
temp4 = sum(temp4,3);
sze = size(im_sub);
beta = betta;
branches = zeros((sze(1)-1)*(sze(2)-1)*2+(sze(1)-1)+(sze(2)-1), 6);
idx = 1;
for y = 1:sze(1)
    for x = 1:sze(2)

        vec_pos = (x-1)*sze(1)+y;
         if y < sze(1)
            vec_pos_d = (x-1)*sze(1)+y+1;
            tempo = gamma*exp(-beta*temp3(y,x));
            branches(idx,:) = [vec_pos,vec_pos_d,0,tempo,tempo,0];
            idx = idx+1;
        end
        
        if x < sze(2)
            vec_pos_r = (x+1-1)*sze(1)+y;
            tempo = gamma*exp(-beta*temp4(y,x));
            branches(idx,:) = [vec_pos,vec_pos_r,0,tempo,tempo,0];
            idx = idx+1;
        end
     end
end

end

end
    
    