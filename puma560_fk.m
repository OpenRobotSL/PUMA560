function [T0t,T061,T062]=puma560_fk(theta1,theta2,theta3,theta4,theta5,theta6)
%广晓模型DH ,theta 2  3 
L6=0; L2=0.2; L7=0; L4=0.248; L5=0.262;L3=0; 
%     theta        d       a-1     alpha-1
MDH=[theta1        0       0        0;
     -theta2-pi/2  0       L6      -pi/2;
     -theta3       L7      L2       0;
     theta4        L4      L3      -pi/2;
     theta5        0       0       -pi/2;
     theta6        0       0        pi/2;
      0            L5      0        pi/2  ];
 
  
  T01=[cos(MDH(1,1))                 -sin(MDH(1,1))               0                 MDH(1,3);
      sin(MDH(1,1))*cos(MDH(1,4))   cos(MDH(1,1))*cos(MDH(1,4))  -sin(MDH(1,4))    -MDH(1,2)*sin(MDH(1,4));
      sin(MDH(1,1))*sin(MDH(1,4))   cos(MDH(1,1))*sin(MDH(1,4))  cos(MDH(1,4))     MDH(1,2)*cos(MDH(1,4));
      0               0                             0                              1];
 
  T12=[cos(MDH(2,1))                 -sin(MDH(2,1))               0                 MDH(2,3);
      sin(MDH(2,1))*cos(MDH(2,4))   cos(MDH(2,1))*cos(MDH(2,4))  -sin(MDH(2,4))    -MDH(2,2)*sin(MDH(2,4));
      sin(MDH(2,1))*sin(MDH(2,4))   cos(MDH(2,1))*sin(MDH(2,4))  cos(MDH(2,4))     MDH(2,2)*cos(MDH(2,4));
      0               0                             0                              1];
 
  T23=[cos(MDH(3,1))                 -sin(MDH(3,1))               0                 MDH(3,3);
      sin(MDH(3,1))*cos(MDH(3,4))   cos(MDH(3,1))*cos(MDH(3,4))  -sin(MDH(3,4))    -MDH(3,2)*sin(MDH(3,4));
      sin(MDH(3,1))*sin(MDH(3,4))   cos(MDH(3,1))*sin(MDH(3,4))  cos(MDH(3,4))     MDH(3,2)*cos(MDH(3,4));
      0               0                             0                              1];
 
  T34=[cos(MDH(4,1))                 -sin(MDH(4,1))               0                 MDH(4,3);
      sin(MDH(4,1))*cos(MDH(4,4))   cos(MDH(4,1))*cos(MDH(4,4))  -sin(MDH(4,4))    -MDH(4,2)*sin(MDH(4,4));
      sin(MDH(4,1))*sin(MDH(4,4))   cos(MDH(4,1))*sin(MDH(4,4))  cos(MDH(4,4))     MDH(4,2)*cos(MDH(4,4));
      0               0                             0                              1];
 
  T45=[cos(MDH(5,1))                 -sin(MDH(5,1))               0                 MDH(5,3);
      sin(MDH(5,1))*cos(MDH(5,4))   cos(MDH(5,1))*cos(MDH(5,4))  -sin(MDH(5,4))    -MDH(5,2)*sin(MDH(5,4));
      sin(MDH(5,1))*sin(MDH(5,4))   cos(MDH(5,1))*sin(MDH(5,4))  cos(MDH(5,4))     MDH(5,2)*cos(MDH(5,4));
      0               0                             0                              1];
 
  T56=[cos(MDH(6,1))                 -sin(MDH(6,1))               0                 MDH(6,3);
      sin(MDH(6,1))*cos(MDH(6,4))   cos(MDH(6,1))*cos(MDH(6,4))  -sin(MDH(6,4))    -MDH(6,2)*sin(MDH(6,4));
      sin(MDH(6,1))*sin(MDH(6,4))   cos(MDH(6,1))*sin(MDH(6,4))  cos(MDH(6,4))     MDH(6,2)*cos(MDH(6,4));
      0               0                             0                              1];

  T6t=[0                 -1              0                 0;
       1                 0               0                 0;
       0                 0               1                 L5;
       0                 0               0                 1];
  
    
 T0t=T01*T12*T23*T34*T45*T56*T6t;
 T061=T01*T12*T23*T34*T45*T56; %逆解需要T06
 T062=T0t*inv(T6t);

end