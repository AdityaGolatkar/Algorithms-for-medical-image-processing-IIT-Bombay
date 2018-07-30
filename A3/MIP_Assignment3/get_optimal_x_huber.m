function [optimal_x,J_min] = get_optimal_x_huber(y,x_n1,x_n2,x_n3,x_n4,alpha,gamma)
tau = 0.005;
iter = 0;
x = (y + x_n1 + x_n2 + x_n3 + x_n4)/5;
tau_min = 1e-7;
while 1
    
    if(iter > 1)
        if(J_curr < J_prev)
            tau = 1.1*tau;
        else
            tau = 0.5*tau;
        end
    end
    
    if(tau < tau_min)
        break;
    end
    
    J_prev = 0;
    
    if(iter > 0)
        J_prev = J_curr;
    end

    J0 = abs(y-x)*abs(y-x);
    D0 = -2*(y-x);

    if(abs(x-x_n1) < gamma)
        J1 = 0.5*abs(x-x_n1)*abs(x-x_n1);
        D1 = (x-x_n1);
    else
        J1 = gamma*abs(x-x_n1) - 0.5*gamma*gamma;
        D1 = gamma*sign(x-x_n1);
    end

    if(abs(x-x_n2) < gamma)
        J2 = 0.5*abs(x-x_n2)*abs(x-x_n2);
        D2 = (x-x_n2);
    else
        J2 = gamma*abs(x-x_n2) - 0.5*gamma*gamma;
        D2 = gamma*sign(x-x_n2);
    end

    if(abs(x-x_n3) < gamma)
        J3 = 0.5*abs(x-x_n3)*abs(x-x_n3);
        D3 = (x-x_n3);
    else
        J3 = gamma*abs(x-x_n3) - 0.5*gamma*gamma;
        D3 = gamma*sign(x-x_n3);
    end

    if(abs(x-x_n4) < gamma)
        J4 = 0.5*abs(x-x_n4)*abs(x-x_n4);
        D4 = (x-x_n4);
    else
        J4 = gamma*abs(x-x_n4) - 0.5*gamma*gamma;
        D4 = gamma*sign(x-x_n4);
    end

    J_curr = (1-alpha)*J0 + alpha*(J1 + J2 + J3 + J4);
    D = (1-alpha)*D0 + alpha*(D1 + D2 + D3 + D4);
    
    %if(iter > 0)
    %    J_diff = abs(J_curr/J_prev - 1);
    %    %if(J_diff < 0.0001)
    %    if(tau < tau_min)
    %        break;
    %    end
    %end
    
    x = x - tau*D;
    iter = iter + 1;
    
end
J_min = J_curr;
optimal_x = x;

