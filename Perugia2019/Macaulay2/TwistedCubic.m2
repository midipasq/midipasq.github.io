--In this file we explore the twisted cubic, which is the image of the 
--map from R^1 to R^3 defined by the parametric equations
--x=t,y=t^2,z=t^3

--Remark: every line that starts with "--" is considered a comment by Macaulay2
--and is not evaluated.  If you are trying to reproduce this code, only enter
--lines which do not start with "--".

--First, define the coordinate ring of affine three space (we'll work
--over the rational numbers since Macaulay2 doesn't work well yet
--over the real or complex numbers)
R=QQ[x,y,z]

--Next, define the coordinate ring of affine one dimensional space
T=QQ[t]

--The parametrization x=t,y=t^2,z=t^3 is realized by the map of coordinate rings
--phi from S to T defined by phi(x)=t, phi(y)=t^2, and phi(z)=t^3.
--We define this map as follows.
phi=map(T,R,{t,t^2,t^3})

--The image of the parametrization is an affine variety in three-dimensional space.
--This means that there are equations that define the image of this map.
--The following command computes the kernel of the ring homorphism phi; this is an
--ideal whose generators cut out the image of the parametrization.
I=ker phi
--Take a moment to check that if you set x=t, y=t^2, and z=t^3, the polynomials 
--generating the ideal will all vanish!

--Intuitively, the twisted cubic is a curve, so it should have dimension one.
--We can recover this with the 'dim' command.
dim I

--Now we explore the tangent variety of the twisted cubic.  This is a surface in
--affine three-space defined by the parametrization x=t+s, y=t^2+2*t*s, z=t^3+3*t^2*s.
--First, define the coordinate ring of affine two-space, which is the 'parameter space.'
S=QQ[s,t]

--Next, define the map and compute its kernel to get the equation defining the tangent surface.
psi=map(S,R,{t+s,t^2+2*t*s,t^3+3*t^2*s})
J=ker psi

--Notice that J is principal; the unique generator is exactly the equation defining the tangent surface!
--As a reality check, find the dimension of J (it should be 2):
dim J

--We have not learned yet how to compute the kernel of a map between rings.
--Computing kernels of maps like this uses 'elimination theory,' which we will
--encounter in chapter 3.