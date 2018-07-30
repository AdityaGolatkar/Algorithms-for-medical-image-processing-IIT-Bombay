function [optimal_x] = get_optimal_x_log(y,S,alpha,gamma)
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
    m1 = abs(p1);
    l1 = log(1+m1/gamma);
    J1 = gamma*sum(m1(:)) - gamma*gamma*sum(l1(:));
    D1 = gamma*sign(p1).*(m1./(m1+gamma));
    
    p2 = x-x_n2;
    m2 = abs(p2);
    l2 = log(1+m2/gamma);
    J2 = gamma*sum(m2(:)) - gamma*gamma*sum(l2(:));
    D2 = gamma*sign(p2).*(m2./(m2+gamma));
    
    p3 = x-x_n3;
    m3 = abs(p3);
    l3 = log(1+m3/gamma);
    J3 = gamma*sum(m3(:)) - gamma*gamma*sum(l3(:));
    D3 = gamma*sign(p3).*(m3./(m3+gamma));
    
    p4 = x-x_n4;
    m4 = abs(p4);
    l4 = log(1+m4/gamma);
    J4 = gamma*sum(m4(:)) - gamma*gamma*sum(l4(:));
    D4 = gamma*sign(p4).*(m4./(m4+gamma));
    
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



