--Get equation for the coupler curve of a 4-Bar linkage.
--lengths of the five bars
r1=2;
r2=1;
r3=3;
r4=4;
r5=4;
--Set up equations.
--Assumes fixed points of (0,0) and (a,0).
a=5;
--other points have variable coordinates given by (x_1,y_1),(x_2,y_2),(s,t)
R=QQ[x_1,x_2,y_1,y_2,s,t]
--define equations giving lengths of bars
f1=x_1^2+y_1^2-r1^2;
f2=(x_2-a)^2+y_2^2-r2^2;
f3=(x_1-x_2)^2+(y_1-y_2)^2-r3^2;
f4=(s-x_1)^2+(t-y_1)^2-r4^2;
f5=(s-x_2)^2+(t-y_2)^2-r5^2;
--define ideal
I=ideal(f1,f2,f3,f4,f5);
--eliminate to get equation in s and t of the coupler curve
J=eliminate(I,{x_1,x_2,y_1,y_2});
--the equation of the coupler curve
F=J_0
--output to a file (otherwise its difficult to copy and paste)
"couplerEq"<<toString(F)<<close