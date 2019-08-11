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
