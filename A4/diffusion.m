function [] = diffusion()
S_0 = 1;
b_0 = 0.1;
g = [1,0.866,0.5,0,-0.5,-0.866;0,0.5,0.866,1,0.866,0.5];
S = [0.5054-0.0217i;0.6874+0.0171i;0.3632+0.1789i;0.3483+0.1385i;0.2606-0.0675i;0.2407+0.1517i];
S = abs(S);
L = [1,0;1,1];
l = ones(3,1);
D = L*L';
tau = 0.01;
lambda = 0.5;


iter = 0;
max_iter = 100;
r = zeros(6,1);
J = zeros(max_iter,1);
d_11 = J;
d_12 = d_11;
d_21 = d_11;
d_22 = d_11;
J_r = zeros(6,3);

while (iter < max_iter)
    
    if(iter > 1)
        if(J_curr < J_prev)
            tau = 1.1*tau;
        else
            tau = 0.5*tau;
        end
    end
    
    if(iter > 0)
        J_prev = J_curr;
    end
    
    
    temp1 = g'*D*g;
    temp2 = b_0*diag(temp1);
    r = S - S_0*exp(-temp2);
    J_curr = r'*r;
    J(iter+1,1) = log(J_curr);
    
    temp3 = g'*L;
    for q=1:6
        J_r(q,1) = S_0*b_0*exp(-temp2(q))*2*temp3(q,1)*g(1,q);
        J_r(q,2) = S_0*b_0*exp(-temp2(q))*2*temp3(q,1)*g(2,q);
        J_r(q,3) = S_0*b_0*exp(-temp2(q))*2*temp3(q,2)*g(2,q);
    end
    H = J_r'*J_r;
    temp4 = diag(H);
    temp5 = diag(temp4);
    
    delta = -tau*(inv(H+lambda*temp5))*J_r'*r; 
    l = l + delta;
    
    L(1,1) = l(1);
    if L(1,1) < 0
        L(1,1) = 0;
    end
    L(1,2) = 0;
    L(2,1) = l(2);
    L(2,2) = l(3);
    if L(2,2) < 0
        L(2,2) = 0;
    end
    
    D = L*L';
    d_11(iter+1,1) = D(1,1);
    d_12(iter+1,1) = D(1,2);
    d_21(iter+1,1) = D(2,1);
    d_22(iter+1,1) = D(2,2);
    
    
    if(iter == 0)
        J_min = J_curr;
    else 
        if(J_curr < J_min)
            J_min = J_curr;
            optimal_D = D;
        end
    end
    
    iter = iter + 1;
end

[eV eD] = eig(D);
disp('D matrix');
D
disp('Principal direction along which the diffusion is strongest');
eV(:,2)
disp('Positive Eigen values to show its PSD matrix')
eD
disp('Strength of diffusion of principal direction as compared to the perpendicular direction');
eD(2,2)/eD(1,1)

figure;
i = 1:max_iter;
plot(i,J);
xlabel('Number Of Iterations');
ylabel('Log of Cost Function');
figure;
plot(i,d_11);
xlabel('Number Of Iterations');
ylabel('D(1,1)');
figure;
plot(i,d_12);
xlabel('Number Of Iterations');
ylabel('D(1,2)');
figure;
plot(i,d_21);
xlabel('Number Of Iterations');
ylabel('D(2,1)');
figure;
plot(i,d_22);
xlabel('Number Of Iterations');
ylabel('D(2,2)');


    
    
    