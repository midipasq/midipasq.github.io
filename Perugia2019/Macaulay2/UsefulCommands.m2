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

--create a function.
--The inputs are in parentheses.
--The output is the last line of the function.
f=(m,n)->(
    k := m+n; 
    k^2
    )
--the last line in your function is what is the output.  Take care not to put
--a semicolon after it!

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

--Testing containment
isSubset(I,J)
isSubset(I,K)
isSubset(K,I)

--(simple) Output to a file.  (It can be useful to output strings to files so that they are formatted nicely for copying and pasting).
R=QQ[x,y,z]
p=x^2+y^2+z^2
"~/Desktop/CouplerEq.txt"<<toString(p)<<close
--If you don't specify the complete path, the file will be saved in the directory you are currenlty working in
"CouplerEq.txt"<<toString(p)<<close

--Loading a file you've created.
--Suppose you've created a file in your local directory called "Perugia.txt"
load("Perugia.txt")