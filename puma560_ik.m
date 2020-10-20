function [thetaF_D,thetaF_R]=puma560_ik_guangxiao(Tbe)
%广晓的DH模型  56有区别
%广晓得逆解，只是针对他自己DH坐标系适用，同时也是考虑了-theta2-90   -theta3，所以在plot，正解时候要注意，取负数
L6=0; L2=0.2; L7=0; L4=0.248; L5=0.262;L3=0; 
%    ti   di     ai-1     alphai-1 
MDH=[0   0       0         0;
     0   0       L6        -pi/2;
     0   L7      L2        0;
     0   L4      L3        -pi/2;
     0   0       0         -pi/2;
     0   0       0         pi/2
     0   L5      0         pi/2];%theta_tool=90  末端工具
px=Tbe(1,4); py=Tbe(2,4); pz=Tbe(3,4);DEG =pi/180;
 rou=sqrt(px^2+py^2);
 
 %theta1

    th11 = atan2(py,px)-atan2(L7/rou,sqrt(1-(L7/rou)^2));
    th12 = atan2(py,px)-atan2(L7/rou,-sqrt(1-(L7/rou)^2));

    
%theta3
    K = (px^2+py^2+pz^2+L6^2-L2^2-L3^2-L4^2-L7^2-2*L6*(px*cos(th11)+py*sin(th11)))/(2*L2);
    %R=sqrt(L3^2+L4^2);

    th31 = atan2(K,sqrt(L3^2+L4^2-K^2))-atan2(L3,L4);       
    th32 = atan2(K,-sqrt(L3^2+L4^2-K^2))-atan2(L3,L4);

    
%theta2 有4个解，分别列出
%th11 th31
Tempy = (L4+L2*sin(th31))*pz-(L3+L2*cos(th31))*(-L6+px*cos(th11)+py*sin(th11));
Tempx = (L3+L2*cos(th31))*pz-(L4+L2*sin(th31))*(L6-px*cos(th11)-py*sin(th11));
thetaT = atan2(Tempy,Tempx);
th21 = thetaT - th31;
%th11 th32
Tempy = (L4+L2*sin(th32))*pz-(L3+L2*cos(th32))*(-L6+px*cos(th11)+py*sin(th11));
Tempx = (L3+L2*cos(th32))*pz-(L4+L2*sin(th32))*(L6-px*cos(th11)-py*sin(th11));
thetaT = atan2(Tempy,Tempx);
th22 = thetaT - th32;
%th12 th31
Tempy = (L4+L2*sin(th31))*pz-(L3+L2*cos(th31))*(-L6+px*cos(th12)+py*sin(th12));
Tempx = (L3+L2*cos(th31))*pz-(L4+L2*sin(th31))*(L6-px*cos(th12)-py*sin(th12));
thetaT = atan2(Tempy,Tempx);
th23 = thetaT - th31;
%th12 th32
Tempy =(L4+L2*sin(th32))*pz-(L3+L2*cos(th32))*(-L6+px*cos(th12)+py*sin(th12));
Tempx = (L3+L2*cos(th32))*pz-(L4+L2*sin(th32))*(L6-px*cos(th12)-py*sin(th12));
thetaT = atan2(Tempy,Tempx);
th24 = thetaT - th32;


%求解theta4，通过th1 th2 th3
ax = Tbe(1,3); ay = Tbe(2,3); az = Tbe(3,3);
%th11 th31 th21
Tempy = ax*sin(th11)-ay*cos(th11);
Tempx = az*cos(th21+th31)-ax*cos(th11)*sin(th31+th21)-ay*sin(th11)*sin(th21+th31);
th41 = atan2(Tempy,Tempx);

%th11 th32 th22
Tempy = ax*sin(th11)-ay*cos(th11);
Tempx = az*cos(th22+th32)-ax*cos(th11)*sin(th32+th22)-ay*sin(th11)*sin(th22+th32);
th42 = atan2(Tempy,Tempx);

%th12 th31 th23
Tempy = ax*sin(th12)-ay*cos(th12);
Tempx = az*cos(th23+th31)-ax*cos(th12)*sin(th31+th23)-ay*sin(th12)*sin(th23+th31);
th43 = atan2(Tempy,Tempx);

%th12 th32 th24
Tempy = ax*sin(th12)-ay*cos(th12);
Tempx = az*cos(th24+th32)-ax*cos(th12)*sin(th32+th24)-ay*sin(th12)*sin(th24+th32);
th44 = atan2(Tempy,Tempx);


%求解theta5 通过th1 th2 th3 th4
%th11 th31 th21 th41
Tempy = ax*(sin(th11)*sin(th41)-sin(th21+th31)*cos(th41)*cos(th11))-ay*(cos(th11)*sin(th41)+sin(th21+th31)*cos(th41)*sin(th11))+az*cos(th21+th31)*cos(th41);
Tempx = az*sin(th21+th31)+ax*cos(th21+th31)*cos(th11)+ay*cos(th21+th31)*sin(th11);
th51 = atan2(Tempy,Tempx);
%th11 th32 th22 th42
Tempy = ax*(sin(th11)*sin(th42)-sin(th22+th32)*cos(th42)*cos(th11))-ay*(cos(th11)*sin(th42)+sin(th22+th32)*cos(th42)*sin(th11))+az*cos(th22+th32)*cos(th42);
Tempx = az*sin(th22+th32)+ax*cos(th22+th32)*cos(th11)+ay*cos(th22+th32)*sin(th11);
th52 = atan2(Tempy,Tempx);
%th12 th31 th23 th43
Tempy = ax*(sin(th12)*sin(th43)-sin(th23+th31)*cos(th43)*cos(th12))-ay*(cos(th12)*sin(th43)+sin(th23+th31)*cos(th43)*sin(th12))+az*cos(th23+th31)*cos(th43);
Tempx = az*sin(th23+th31)+ax*cos(th23+th31)*cos(th12)+ay*cos(th23+th31)*sin(th12);
th53 = atan2(Tempy,Tempx);
%th12 th32 th24 th44
Tempy = ax*(sin(th12)*sin(th44)-sin(th24+th32)*cos(th44)*cos(th12))-ay*(cos(th12)*sin(th44)+sin(th24+th32)*cos(th44)*sin(th12))+az*cos(th24+th32)*cos(th44);
Tempx = az*sin(th24+th32)+ax*cos(th24+th32)*cos(th12)+ay*cos(th24+th32)*sin(th12);
th54 = atan2(Tempy,Tempx);

%求解theta6
nx = Tbe(1,1); ny = Tbe(2,1); nz = Tbe(3,1);
%th11 th31 th21 th41 th51
Tempy = nx*(cos(th41)*sin(th11)+sin(th21+th31)*sin(th41)*cos(th11))-ny*(cos(th11)*cos(th41)-sin(th21+th31)*sin(th11)*sin(th41))-nz*cos(th21+th31)*sin(th41);
Tempx1 = -nx*(cos(th21+th31)*cos(th11)*sin(th51)-cos(th51)*sin(th11)*sin(th41)+sin(th21+th31)*cos(th11)*cos(th41)*cos(th51));
Tempx2 = -ny*(cos(th21+th31)*sin(th11)*sin(th51)+cos(th11)*cos(th51)*sin(th41)+sin(th21+th31)*cos(th41)*cos(th51)*sin(th11));
Tempx3 = -nz*(sin(th21+th31)*sin(th51)-cos(th21+th31)*cos(th41)*cos(th51));
Tempx=Tempx1+Tempx2+Tempx3;
th61 = atan2(Tempy,Tempx);
%th11 th32 th22 th42 th52
Tempy = nx*(cos(th42)*sin(th11)+sin(th22+th32)*sin(th42)*cos(th11))-ny*(cos(th11)*cos(th42)-sin(th22+th32)*sin(th11)*sin(th42))-nz*cos(th22+th32)*sin(th42);
Tempx1 = -nx*(cos(th22+th32)*cos(th11)*sin(th52)-cos(th52)*sin(th11)*sin(th42)+sin(th22+th32)*cos(th11)*cos(th42)*cos(th52));
Tempx2 = -ny*(cos(th22+th32)*sin(th11)*sin(th52)+cos(th11)*cos(th52)*sin(th42)+sin(th22+th32)*cos(th42)*cos(th52)*sin(th11));
Tempx3 = -nz*(sin(th22+th32)*sin(th52)-cos(th22+th32)*cos(th42)*cos(th52));
Tempx = Tempx1 + Tempx2 + Tempx3;
th62 = atan2(Tempy,Tempx);
%th12 th31 th23 th43 th53
Tempy = nx*(cos(th43)*sin(th12)+sin(th23+th31)*sin(th43)*cos(th12))-ny*(cos(th12)*cos(th43)-sin(th23+th31)*sin(th12)*sin(th43))-nz*cos(th23+th31)*sin(th43);
Tempx1 = -nx*(cos(th23+th31)*cos(th12)*sin(th53)-cos(th53)*sin(th12)*sin(th43)+sin(th23+th31)*cos(th12)*cos(th43)*cos(th53));
Tempx2 = -ny*(cos(th23+th31)*sin(th12)*sin(th53)+cos(th12)*cos(th53)*sin(th43)+sin(th23+th31)*cos(th43)*cos(th53)*sin(th12));
Tempx3 = -nz*(sin(th23+th31)*sin(th53)-cos(th23+th31)*cos(th43)*cos(th53));
Tempx=Tempx1+Tempx2+Tempx3;
th63 = atan2(Tempy,Tempx);
%th12 th32 th24 th44 th54
Tempy = nx*(cos(th44)*sin(th12)+sin(th24+th32)*sin(th44)*cos(th12))-ny*(cos(th12)*cos(th44)-sin(th24+th32)*sin(th12)*sin(th44))-nz*cos(th24+th32)*sin(th44);
Tempx1 = -nx*(cos(th24+th32)*cos(th12)*sin(th54)-cos(th54)*sin(th12)*sin(th44)+sin(th24+th32)*cos(th12)*cos(th44)*cos(th54));
Tempx2 = -ny*(cos(th24+th32)*sin(th12)*sin(th54)+cos(th12)*cos(th54)*sin(th44)+sin(th24+th32)*cos(th44)*cos(th54)*sin(th12));
Tempx3 = -nz*(sin(th21+th32)*sin(th54)-cos(th24+th32)*cos(th44)*cos(th54));
Tempx=Tempx1+Tempx2+Tempx3;
th64 = atan2(Tempy,Tempx);

%最终结果，矩阵显示
thetaF = zeros(8,6);
thetaF(1,1) = th11/DEG; thetaF(1,2) = th21/DEG; thetaF(1,3) = th31/DEG; thetaF(1,4) = th41/DEG; thetaF(1,5) = th51/DEG; thetaF(1,6) = th61/DEG;
thetaF(2,1) = th11/DEG; thetaF(2,2) = th22/DEG; thetaF(2,3) = th32/DEG; thetaF(2,4) = th42/DEG; thetaF(2,5) = th52/DEG; thetaF(2,6) = th62/DEG;
thetaF(3,1) = th12/DEG; thetaF(3,2) = th23/DEG; thetaF(3,3) = th31/DEG; thetaF(3,4) = th43/DEG; thetaF(3,5) = th53/DEG; thetaF(3,6) = th63/DEG;
thetaF(4,1) = th12/DEG; thetaF(4,2) = th24/DEG; thetaF(4,3) = th32/DEG; thetaF(4,4) = th44/DEG; thetaF(4,5) = th54/DEG; thetaF(4,6) = th64/DEG;
thetaF(5,1) = th11/DEG; thetaF(5,2) = th21/DEG; thetaF(5,3) = th31/DEG; thetaF(5,4) = th41/DEG+180; thetaF(5,5) = -th51/DEG; thetaF(5,6) = th61/DEG+180;
thetaF(6,1) = th11/DEG; thetaF(6,2) = th22/DEG; thetaF(6,3) = th32/DEG; thetaF(6,4) = th42/DEG+180; thetaF(6,5) = -th52/DEG; thetaF(6,6) = th62/DEG+180;
thetaF(7,1) = th12/DEG; thetaF(7,2) = th23/DEG; thetaF(7,3) = th31/DEG; thetaF(7,4) = th43/DEG+180; thetaF(7,5) = -th53/DEG; thetaF(7,6) = th63/DEG+180;
thetaF(8,1) = th12/DEG; thetaF(8,2) = th24/DEG; thetaF(8,3) = th32/DEG; thetaF(8,4) = th44/DEG+180; thetaF(8,5) = -th54/DEG; thetaF(8,6) = th64/DEG+180;



%为了画图,maltab画图里面q2  q3是正的，逆解出来也是正的，但是DH建模是负，所以画图需要负的，如果不画图用光效整洁，就需要再次取反
thetaF(:,2)=-thetaF(:,2);
thetaF(:,3)=-thetaF(:,3);

thetaF_R=thetaF/180*pi;

thetaF_D=thetaF;

end