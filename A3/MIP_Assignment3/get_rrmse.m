function [rrmse] = get_rrmse(img1,img2)
Ap = abs(img1);
Bp = abs(img2);
Ap = Ap(:);
Bp = Bp(:);
Cp = Ap - Bp;
Cp2 = Cp.*Cp;
Ap2 = Ap.*Ap;
S_Cp2 = sum(Cp2);
S_Ap2 = sum(Ap2);
rrmse = sqrt(S_Cp2/S_Ap2);

