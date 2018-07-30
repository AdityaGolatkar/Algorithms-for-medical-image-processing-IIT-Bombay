function [branches] = edges(inp_crop,gamma)

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
beta = 1/((2*temp5)/(prod(sztest(1:2)))-sztest(1)-sztest(2));

branches = zeros((sztest(1)-1)*(sztest(2)-1)*2+(sztest(1)-1)+(sztest(2)-1),6);

%temp6 = double(gamma*exp(-beta*(sum(temp3,3))));
%temp7 = double(gamma*exp(-beta*(sum(temp4,3))));
diff1 = zeros(1,3);
ct = 1;
for y = 1:sz(1)
    for x = 1:sz(2)
        pos = (x-1)*sz(1)+y;
        now = get_rgb_double(inp_crop,x,y);
        if x<sz(2)
            curr = double(inp_crop(y,x,:));
            nex  = double(inp_crop(y,x+1,:));
            diff = double(curr-nex);
            diff1(1,1) = diff(1,1,1);
            diff1(1,2) = diff(1,1,2);
            diff1(1,3) = diff(1,1,3);
            next = get_rgb_double(inp_crop,x+1,y);
            %temp6 = gamma*exp(-beta*(norm(now-next)^2));
            %temp6 = norm(now-next)^2;
            temp6 = sum(now);
            branches(ct,:)=[pos,pos+sz(1),0,temp6,temp6,0];
           
            ct=ct+1;
        end
        if y<sz(1)
            curr = double(inp_crop(y,x,:));
            nex  = double(inp_crop(y+1,x,:));
            diff = double(curr-nex);
            diff1(1,1) = diff(1,1,1);
            diff1(1,2) = diff(1,1,2);
            diff1(1,3) = diff(1,1,3);
            next = get_rgb_double(inp_crop,x,y+1);
            %temp6 = gamma*exp(-beta*(norm(now-next)^2));
            %temp6 = norm(now-next)^2;
            branches(ct,:)=[pos,pos+1,0,temp6,temp6,0];
           
            ct=ct+1;
        end
    end
end