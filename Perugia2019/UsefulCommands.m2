--Useful commands to know in Macaulay2:

--Lists
L={3,4,5,9}
--get the length of the list
length L
--or
#L
--Access 1st element of list
L_0
--or
L#0
--or
first L
--Access last element of list
L_3
--or
L#3
--or
last L

--Nested lists
L={{3,4,5},{6,5,2}}
--Access the '4'
L_0_1
--Access the '6'
L_1_0

--creating a matrix
M=matrix({{3,4,5},{6,5,2}})
--access the '3'
M_(0,0)
--access the '6'
M_(1,0)
--multiply matrices
N=matrix({{2,3},{4,2},{6,7}})
M*N
--take transpose
transpose M
--get the list of rows of M
entries M

--Make a matrix into a list
M=matrix({{3,4,5},{6,5,2}})
E=entries M --returns the list {{3,4,5},{6,5,2}}
flatten E --returns the list {3,4,5,6,5,2} (this removes one set of inner braces)

--create a function.
--The inputs are in parentheses.
--The output is the last line of the function.
f=(m,n)->(
    k := m+n; 
    k^2
    )
--the last line in your function is what is the output.  Take care not to put
--a semicolon after it!

--Apply a function to a list
L={0,1,2,3}
apply(L,i->i^2)--squares every element of the list

--Find a Grobner basis
R=QQ[x,y,MonomialOrder=>Lex] --Default is GRevLex.  Other standard option is GLex.
I=ideal(x^2+y^2-1,x^3+y^2-1)
gens gb I

--Specify a monomial order by weights
R=QQ[x,y,MonomialOrder=>{Weights=>{1,0},Weights=>{1,0}}]--defines Lex order
--get all monomials of degree 3
L=flatten entries basis(3,R)
--sort them according to your chosen monomial order
sort L

--Testing equality
R=QQ[x,y]
I=ideal(x^2,y^2)
J=ideal(x^2,y^2-x^2)
K=ideal(x^2,x^2-y)
I==K

--Testing containment (using ideals defined above)
isSubset(I,J)
isSubset(I,K)
isSubset(K,I)

--get the generators of an ideal
R=QQ[x,y]
I=ideal(x^2,y^2)
gens I --returns the matrix matrix {{x^2,y^2}}
flatten entries gens I --returns the list {x^2,y^2}

--Intersect, quotient, product, sum, and saturate
R=QQ[x,y]
m1=ideal(x,y)
m2=ideal(x-1,y)
I=intersect(m1,m2) --intersects the ideals m1 and m2
J=intersect(m1^2,m2^2) --intersects the squares of the ideals m1 and m2
K=m1*m2 --takes product of the ideals m1 and m2
J:m1 --takes the colon J:m1
J:I --takes the colon J:I
saturate(J,m1) --computes the saturation J:m1^infinity
K=m1+m2 --takes the sum of m1 and m2

--Get a monomial basis for the quotient of a zero-dimensional ideal
S=QQ[x,y]
I=ideal(x^3*y-x,y^2-5*y*x+6)
basis(S/I)--gives all monomials outside LT(I) with respect to monomial order
--Try a different monomial order
S=QQ[x,y,MonomialOrder=>Lex]
I=ideal(x^3*y-x,y^2-5*y*x+6)
basis(S/I)--notice the basis has changed due to the change in monomial order

--Get a monomial basis in a fixed degree for the quotient by any ideal
S=QQ[x,y,z]
I=ideal(x^2-y,x^3-z)
basis(4,S/I)--gives monomials of degree 4 outside of LT(I) with respect to monomial order

--Output a string to a file.  (It can be useful to output strings to files so that they are formatted nicely for copying and pasting).
R=QQ[x,y,z]
p=x^2+y^2+z^2
"~/Desktop/CouplerEq.txt"<<toString(p)<<close
--If you don't specify the complete path, the file will be saved in the directory you are currenlty working in
"CouplerEq.txt"<<toString(p)<<close

--Loading a file you've created.
--Suppose you've created a file in your local directory called "Perugia.txt"
load("Perugia.txt")


