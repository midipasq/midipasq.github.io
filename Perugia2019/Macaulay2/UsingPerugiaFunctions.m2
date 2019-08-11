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


--Using the multiplicationMatrix command
S=QQ[x,y]
I=intersect(ideal(x,y),ideal(x-1,y),ideal(x,y-1))--define an ideal
f=x^2+y^2--define a polynomial to multiply by
M=multiplicationMatrix(f,I)
eigenvalues(M) --gives the values of f evaluated at the points

--Using the command divisionAlgorithm
