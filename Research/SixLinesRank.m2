--Polynomials f1,f2,f3,f4,f5, and f6 are the polynomials in Theorem 5.1*
--The code below produces 
--Verify that the generators of the apolar ideal in degree three and less define 
--a zero-dimensional scheme (for f1,f2,f3,f4,f5 it is a variety of six points).
--For f1,f2,f3,f4, and f5, the points can be identified using the decompose command.

--f1
--define field containing coefficients of linear factors in f1
F=toField(ZZ[t]/ideal(t^2+t+1))
R=F[x,y,z]
f1=x*y*z*(x+y+z)*(x+t*y+t^2*z)*(x+t^2*y+t*z)
--the product is actually defined over the rationals, so now we switch rings
S=QQ[x,y,z]
f1=sub(f1,S)
--the command inverseSystem produces the apolar ideal
I1=inverseSystem f1
--select generators of I1 of degree three and less
P1=ideal select(I1_*,f->degree f=={3})
(radical P1)==P1--P1 is a radical ideal
decompose P1--decomposes into six points, which gives Waring decomposition.

--f2
S=QQ[x,y,z]
f2=x*y*z*(x+y+z)*(x^2 + x*y + 1/3*y^2 + x*z + y*z + z^2)
I2=inverseSystem f2
P2=ideal select(I2_*,f->(first degree f)<4)
P2==radical P2
decompose P2--decomposes into six points, which give Waring decomposition

--f3 (f3 is not defined over the rationals, so we need to introduce an extra variable)
--E is a root of x^2-x+1/3
S=QQ[E,x,y,z]
--f3=x*y*z*(x+y+z)*(x+baro*y+o*z)*(x+barE*y+E*z) o is a root of x^2 - x + 1, baro is its conjugate, barE is conjugate of E
--minimal polynomial of E
MP=E^2-E+1/3;
--Annihilating cubics (produced using Sage)
H1=(-6*E+2)*x*y^2+(6*E-3)*y^3-(6*E-3)*z^3+x^3+2*x*y*z+(6*E-4)*x*z^2+(3*E-3)*y*z^2;
H2=(3*E-1)*y^3-3*E*x*y^2+x^2*y;
H3=(-3*E+2)*z^3+(3*E-3)*x*z^2+x^2*z;
H4=y^2*z+(3*E-2)*y*z^2;
--Apolar algebra of f5 in degree 3 (we throw in the minimal polynomial of E since it also vanishes)
I3=ideal(MP,H1,H2,H3,H4)
decompose I3--this consists of 6 ideals -- a little inspection reveals these are the points in Table 2 of the paper
--We produce a new ring to verify that H1,H2,H3,H4 actually annihilate f3
K=toField(QQ[E]/ideal(E^2-E+1/3))
R=K[x,y,z]
H1=(-6*E+2)*x*y^2+(6*E-3)*y^3-(6*E-3)*z^3+x^3+2*x*y*z+(6*E-4)*x*z^2+(3*E-3)*y*z^2;
H2=(3*E-1)*y^3-3*E*x*y^2+x^2*y;
H3=(-3*E+2)*z^3+(3*E-3)*x*z^2+x^2*z;
H4=y^2*z+(3*E-2)*y*z^2;
--Use the fact that barE=1-E, o=3E-1, and baro=2-3E to define f3 in the ring R
f3=x*y*z*(x+y+z)*(x+(2-3*E)*y+(3*E-1)*z)*(x+(1-E)*y+E*z)
--The following all return 0, so H1,H2,H3,H4 all annihilate f3
diff(H1,f3)
diff(H2,f3)
diff(H3,f3)
diff(H4,f3)

--f4
--E is a root of x^2-x+1/3
S=QQ[E,x,y,z]
--f4=x*y*z*(x+y+z)*(x+o*y+baro*z)*(x+E*y+barE*z) o is a root of x^2 - x + 1, baro is its conjugate, barE is conjugate of E
--This is conjugate of the third, so has the same behavior

--f5 (A3 arrangement)
S=QQ[x,y,z]
f5=x*y*z*(x+y+z)*(x+y)*(y+z)
I5=inverseSystem f5
P5=ideal select(I5_*,f->(first degree f)<4)
--verify P5 defines a variety
P5==radical P5
decompose P5--decomposes into six points, which gives Waring decomposition.

--f6
S=QQ[x,y,z]
f6=x^3*y*z*(x+y+z)
I6=inverseSystem f6
P6=ideal select(I6_*,f->(first degree f)<4)
--verify P6 is a saturated ideal
P6==saturate(P6)
--verify P6 defines a 0-dimensional scheme
codim P6
