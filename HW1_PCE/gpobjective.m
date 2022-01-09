function f = gpobjective(x,Y)
xi=[0 0 0 0 ;
0 0 1.732 0 ;
0 0 0 1.732 ;
0 0 -1.732 0 ;
0 0 0 -1.732 ;
0 0 1.732 -1.732 ;
0 0 -1.732 1.732 ;
0 0 1.732 1.732 ;
0 0 -1.732 -1.732 ;
1.732 0 0 0 ;
1.732 0 1.732 0 ;
1.732 0 0 1.732 ;
1.732 0 -1.732 0 ;
1.732 0 0 -1.732 ;
1.732 0 1.732 -1.732 ;
1.732 0 -1.732 1.732 ;
1.732 0 1.732 1.732 ;
1.732 0 -1.732 -1.732 ;
0 1.732 0 0 ;
0 1.732 1.732 0 ;
0 1.732 0 1.732 ;
0 1.732 -1.732 0 ;
0 1.732 0 -1.732 ;
0 1.732 1.732 -1.732 ;
0 1.732 -1.732 1.732 ;
0 1.732 1.732 1.732 ;
0 1.732 -1.732 -1.732 ;
-1.732 0 0 0 ;
-1.732 0 1.732 0 ;
-1.732 0 0 1.732 ;
-1.732 0 -1.732 0 ;
-1.732 0 0 -1.732 ;
-1.732 0 1.732 -1.732 ;
-1.732 0 -1.732 1.732 ;
-1.732 0 1.732 1.732 ;
-1.732 0 -1.732 -1.732 ;
0 -1.732 0 0 ;
0 -1.732 1.732 0 ;
0 -1.732 0 1.732 ;
0 -1.732 -1.732 0 ;
0 -1.732 0 -1.732 ;
0 -1.732 1.732 -1.732 ;
0 -1.732 -1.732 1.732 ;
0 -1.732 1.732 1.732 ;
0 -1.732 -1.732 -1.732 ;
1.732 -1.732 0 0 ;
1.732 -1.732 1.732 0 ;
1.732 -1.732 0 1.732 ;
1.732 -1.732 -1.732 0 ;
1.732 -1.732 0 -1.732 ;
1.732 -1.732 1.732 -1.732 ;
1.732 -1.732 -1.732 1.732 ;
1.732 -1.732 1.732 1.732 ;
1.732 -1.732 -1.732 -1.732 ;
-1.732 1.732 0 0 ;
-1.732 1.732 1.732 0 ;
-1.732 1.732 0 1.732 ;
-1.732 1.732 -1.732 0 ;
-1.732 1.732 0 -1.732 ;
-1.732 1.732 1.732 -1.732 ;
-1.732 1.732 -1.732 1.732 ;
-1.732 1.732 1.732 1.732 ;
-1.732 1.732 -1.732 -1.732 ;
1.732 1.732 0 0 ;
1.732 1.732 1.732 0 ;
1.732 1.732 0 1.732 ;
1.732 1.732 -1.732 0 ;
1.732 1.732 0 -1.732 ;
1.732 1.732 1.732 -1.732 ;
1.732 1.732 -1.732 1.732 ;
1.732 1.732 1.732 1.732 ;
1.732 1.732 -1.732 -1.732 ;
-1.732 -1.732 0 0 ;
-1.732 -1.732 1.732 0 ;
-1.732 -1.732 0 1.732 ;
-1.732 -1.732 -1.732 0 ;
-1.732 -1.732 0 -1.732 ;
-1.732 -1.732 1.732 -1.732 ;
-1.732 -1.732 -1.732 1.732 ;
-1.732 -1.732 1.732 1.732 ;
-1.732 -1.732 -1.732 -1.732 ;
];


for i=1:81
    for j=1:81
        k(i,j)= x(5)^2*exp(-0.5*sqrt((xi(i,1)-xi(j,1))^2/x(1)^2+(xi(i,2)-xi(j,2))^2/x(2)^2+...
        (xi(i,3)-xi(j,3))^2/x(3)^2+(xi(i,4)-xi(j,4))^2/x(4)^2));
    end
end
f=-0.5*Y'*inv(k)*Y-0.5*log(det(k))-81/2*log(2*3.14);
% argmax_x f = argmin_x f
f=-f;