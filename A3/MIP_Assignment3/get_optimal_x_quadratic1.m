function [optimal_x] = get_optimal_x_quadratic(y,x_n1,x_n2,x_n3,x_n4,alpha)
% x_n1:upper neighbor, x_n2:lower neighbor, x_n3:left neighbor, x_n4:right neighbor
optimal_x = (1-alpha)*y + alpha*(x_n1+x_n2+x_n3+x_n4);
optimal_x = optimal_x/(1+3*alpha);
