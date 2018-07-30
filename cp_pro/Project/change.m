function [] = change(img)
c1 = img(:,:,1);
c2 = img(:,:,2);
c3 = img(:,:,3);
f1 = find(c1 <5);
f2 = find(c2 <5);
f3 = find(c3 <5);
cm = intersect(f1,f2);
cm = intersect(cm,f3);
c1(cm) = 0;
c2(cm) = 255;
c3(cm) = 0;
res(:,:,1) = c1;
res(:,:,2) = c2;
res(:,:,3) = c3;
imshow(res);
imwrite(res,'h1.jpg');