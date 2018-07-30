function [optimal_x] = get_optimal_x_huber(y,S,alpha,gamma)
x = ifft2(y);
tau = 0.005;
iter = 0;
iter_max = 200;
J = zeros(iter_max,1);
J_min = 0;
while (iter < iter_max) 
    
    if(iter > 1)
        if(J_curr < J_prev)
            tau = 1.1*tau;
        else
            tau = 0.5*tau;
        end
    end
    
    J_prev = 0;
    
    if(iter > 0)
        J_prev = J_curr;
    end

    m0 = abs(y-(S.*fft2(x))).^2;
    J0 = sum(m0(:));
    D0 = -2*ifft2(S'*(y-(S.*fft2(x))));

    x_n1 = circshift(x,[1 0]);
    x_n2 = circshift(x,[-1 0]);
    x_n3 = circshift(x,[0 1]);
    x_n4 = circshift(x,[0 -1]);
    
    m1 = abs(x-x_n1);
    p1 = m1;
    m2 = abs(x-x_n2);
    p2 = m2;
    m3 = abs(x-x_n3);
    p3 = m3;
    m4 = abs(x-x_n4);
    p4 = m4;
    D1 = zeros(size(x));
    D2 = zeros(size(x));
    D3 = zeros(size(x));
    D4 = zeros(size(x));
    
    b1 = (m1 <= gamma);
    D1(b1) = x(b1)-x_n1(b1);
    D1(~b1) = gamma*sign(x(~b1)-x_n1(~b1));
    p1(b1) = 0.5*m1(b1).^2;
    p1(~b1) = gamma*m1(~b1)-0.5*gamma*gamma;
    J1 = sum(p1(:));
     
    b2 = (m2 <= gamma);
    D2(b2) = x(b2)-x_n2(b2);
    D2(~b2) = gamma*sign(x(~b2)-x_n2(~b2));
    p2(b2) = 0.5*m2(b2).^2;
    p2(~b2) = gamma*m2(~b2)-0.5*gamma*gamma;
    J2 = sum(p2(:));
    
    b3 = (m3 <= gamma);
    D3(b3) = x(b3)-x_n3(b3);
    D3(~b3) = gamma*sign(x(~b3)-x_n3(~b3));
    p3(b3) = 0.5*m3(b3).^2;
    p3(~b3) = gamma*m3(~b3)-0.5*gamma*gamma;
    J3 = sum(p3(:));
    
    b4 = (m4 <= gamma);
    D4(b4) = x(b4)-x_n4(b4);
    D4(~b4) = gamma*sign(x(~b4)-x_n4(~b4));
    p4(b4) = 0.5*m1(b4).^2;
    p4(~b4) = gamma*m4(~b4)-0.5*gamma*gamma;
    J4 = sum(p4(:));
    
    J_curr = (1-alpha)*J0 + alpha*(J1 + J2 + J3 + J4);
    D = (1-alpha)*D0 + alpha*(D1 + D2 + D3 + D4);    
    
    x = x - tau*D;
    J(iter+1,1) = J_curr;
    
    if(iter == 0)
        J_min = J_curr;
    else 
        if(J_curr < J_min)
            J_min = J_curr;
            optimal_x = x;
        end
    end
    
    iter = iter + 1;
    
end

%optimal_x = x;
%iter
iterations = 1:iter;
figure;
plot(iterations,J(iterations),'-o');
xlabel('Number Of Iterations');
ylabel('Negative Log-Likelihood');

%figure;
%imshow(optimal_x);
