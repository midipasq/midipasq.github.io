--Here's a short explanation about using the functions in the file "PerugiaFunctions.txt"

--First, download the file and save it somewhere.
--I've saved the file at ~\Desktop\PerugiaFunctions.txt
--The following command loads the file.
load "~/Desktop/PerugiaFunctions.txt"

--Using the rref command (over QQ)
M=matrix(QQ,{{2,3,4},{3,2,5}})--define a matrix (make sure its over a field)
rref M

--Using the rref command (over a finite field)
A=(ZZ/3)[x]/ideal(x^3+x^2+2)--make sure the polynomial is irreducible
F= GF A --tells Macaulay2 this is a finite field
M=matrix(F,{{1, x^2, 1},{0,x,2}})
rref M

--Using the rref command (over a rational function field)
F=frac(QQ[x,y])
M=matrix(F,{{1,x,x^2},{x+1,x^2+2,x+y}})
rref M

--Using the rref command (over fraction field of a domain)
F=frac(QQ[x,y]/ideal(y^2-x^3-x-1))
F=frac(ZZ/7[x,y]/ideal(y^2-x^3-x-1))
M=matrix(F,{{1,x,x^2},{x+1,x^2+2,x+y}})
rref M

--Using the multiplicationMatrix command
S=QQ[x,y]
I=intersect(ideal(x,y),ideal(x-1,y),ideal(x,y-1))--define an ideal
f=x^2+y^2--define a polynomial to multiply by
(B,M)=multiplicationMatrix(f,I)
B1=basis(S/I);
f*B1--compare columns of multiplication matrix to entries of this
eigenvalues(M) --gives the values of f evaluated at the points of V(I)

--Using the command divisionAlgorithm
S=QQ[x,y,z,MonomialOrder=>GLex]
f=x^3-x^2*y-x^2*z+x;
f1=x^2*y-z;
f2=x*y-1;
L={f1,f2};
(Q1,R1)=divisionAlgorithm(f,L)
(Q2,R2)=divisionAlgorithm(f,reverse L)--different orders can have different remainders

--Using the command sPolynomial
S=QQ[x,y,z,MonomialOrder=>GLex]
f1=x^2*y-z;
f2=x*y-1;
sPolynomial(f1,f2)--notice this is the same as R1-R2

--Using the command buchbergerAlgorithm
S=QQ[x,y,z,MonomialOrder=>GLex]
f1=x^2*y-z;
f2=x*y-1;
L={f1,f2};
buchbergerAlgorithm(L)--notice this is highly redundant
flatten entries gens gb ideal L--compare to this output

--Using couplerCurve
R=QQ[x,y]
a=5; 
r1=1; 
r2=2;
r3=4;
r4=5;
r5=4;
couplerCurve(a,r1,r2,r3,r4,r5,R)

--Using tringImplicitization
R=QQ[x,y]
n=1;
d=2;
f=trigImplicitization(n,d,R)

--Using getRealPlanarSingularPointsNAG
R=QQ[x,y]
n=1;
d=2;
f=trigImplicitization(n,d,R)
I=ideal(x,y)--we know the origin is a highly singular point
getRealPlanarSingularPointsNAG(f,I)

--Using getTrigSingularPointsNAG
n=1;
d=2;
filename="~/Desktop/Trig.txt"
getTrigSingularPointsNAG(n,d,filename)

--Using euclideanDistanceDegree
S=QQ[x,y]
f=(x^2+y^2+x)^2-(x^2+y^2);--equation for a cardioid
L={{-1,0},{2,3},{0,1},{0,2},{3,4}}--list of points to compute distance from
scan(L,P->print euclideanDistanceDegree(ideal f,P))--some give a value of 2
R=entries random(QQ^10,QQ^2)
scan(R,P->print euclideanDistanceDegree(ideal f,P))--probably all give a value of 3
euclideanDistanceDegree(ideal f,P)

f=trigImplicitization(1,2,S)--trig curve
scan(R,P->print euclideanDistanceDegree(ideal f,P))--probably all give 6
scan(L,P->print euclideanDistanceDegree(ideal f,P))--two give 4
