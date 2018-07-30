function [branches] = edges123(inp_crop,gamma)

sz = size(inp_crop);

temp1 = circshift(inp_crop,[-1,0]);
temp2 = circshift(inp_crop,[0,-1]);

temp3 = bsxfun(@minus,temp1,inp_crop);
temp3 = temp3.*temp3;
temp4 = bsxfun(@minus,temp2,inp_crop);
temp4 = temp4.*temp4;

temp5 = sum(sum(sum(temp3)))+sum(sum(sum(temp4)));
beta = 1/((2*temp5)/prod(sz(1:2)));

branches = zeros((sz(1)-1)*(sz(2)-1)*2+(sz(1)-1)+(sz(2)-1),6);

temp6 = double(gamma*exp(-beta*(sum(temp3,3))));
temp7 = double(gamma*exp(-beta*(sum(temp4,3))));

ct = 1;
for y = 1:sz(1)
    for x = 1:sz(2)
        pos = (x-1)*sz(1)+y;
        
        if x<sz(2)
            
            branches(ct,:)=[pos,pos+sz(1),0,temp6(y,x),temp6(y,x),0];
           
            ct=ct+1;
        end
        if y<sz(1)
            
            branches(ct,:)=[pos,pos+1,0,temp7(y,x),temp7(y,x),0];
           
            ct=ct+1;
        end
    end
end

        
        
        
