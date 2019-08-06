--Implement integer programming using Grobner bases.


--Cones
R=QQ[x,y,w_1,w_2,w_3,MonomialOrder=>{Weights=>{1,1,0,0,0},Weights=>{0,0,1,22,15}}]
I=ideal(w_1-x^2*y^5,w_2-x^3*y,w_3-x^3*y^3)
m=x^100*y^100
m % I

R=QQ[w_1,w_2,w_3,Degrees=>{{2,5},{3,1},{3,3}}]
basis({100,100},R)
