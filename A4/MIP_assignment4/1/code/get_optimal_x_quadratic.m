function [optimal_x] = get_optimal_x_quadratic(y,S,alpha)
x = ifft2(y);
tau = 0.005;
iter = 0;
max_iter = 200;
J = zeros(max_iter,1);
J_min = 0;
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
    
    p0 = y-(S.*fft2(x));
    m0 = abs(p0).^2;
    J0 = sum(m0(:));
    D0 = -2*ifft2((S'*p0));
    %D0 = -2*ifft2(S'.*p0);
    
    x_n1 = circshift(x,[1 0]);
    x_n2 = circshift(x,[-1 0]);
    x_n3 = circshift(x,[0 1]);
    x_n4 = circshift(x,[0 -1]);
    
    p1 = x-x_n1;
    m1 = abs(p1).^2;
    J1 = sum(m1(:));
    D1 = 2*p1;
    
    p2 = x-x_n2;
    m2 = abs(p2).^2;
    J2 = sum(m2(:));
    D2 = 2*p2;
    
    p3 = x-x_n3;
    m3 = abs(p3).^2;
    J3 = sum(m3(:));
    D3 = 2*p3;
    
    p4 = x-x_n4;
    m4 = abs(p4).^2;
    J4 = sum(m4(:));
    D4 = 2*p4;
    
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


