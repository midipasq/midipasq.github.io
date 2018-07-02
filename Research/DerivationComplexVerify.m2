--File for verifying examples in the paper "Homological Obstructions to Freeness of Multi-Arrangements" --
--package "HyperplaneArrangements" is automatically loaded in the following file
load "DerivationComplex.m2" --change to path where file is located, for instance load "~/Desktop/M2/DerivationComplex.m2"

-----------------------------------------------------------------
----Example 1.2: Zariski closed set of free multiplicities ------
-----------------------------------------------------------------
S=QQ[x,y,z];
--set up the arrangement in Example 1.2
a1=3; a2=-3;
L={x,y,z,x-a1*z,x-a2*z,y-z}; A=arrangement L; m={3,3,3,1,1,3};
--(A,m) is free if and only if a1/a2 is a root of unity (over reals,
--if and only if a1=-a2
prune image der(A,m)

-----------------------------------------------------------------
----Example 1.3: Totally non-free is not combinatorial ----------
-----------------------------------------------------------------
S=QQ[x,y,z];
--set up the arrangement in Example 1.3
--be sure neither a1 or a2 is equal to 1; this will introduce another triple point
a1=3; a2=-3;
L={x,y,z,x-a1*y,x-a2*y,y-z,z-x}; A=arrangement L; m={3,3,3,1,1,1,1};
--if a1=-a2 and neither is equal to 1 then D(A,m) is free, otherwise not--
prune image der(A,m)

--for slightly more interesting example, set up over a finite field (result is not guaranteed here, but holds in many cases)--
F=ZZ/7; S=F[x,y,z];
--again, be sure a1 is not equal to 1
a1=3; a2=5;
L={x,y,z,x-a1*y,x-a2*y,y-z,z-x}; A=arrangement L; 
--a1^3=a2^3=-1, which is not 1, so the following multiplicity will be free
m={4,4,4,1,1,1,1};
prune image der(A,m)

---------------------------------------------------------------------------
----Example 1.4: Further counter-examples to Orlik's conjecture -----------
---------------------------------------------------------------------------
--create arrangement of rank r+1 in Example 1.4--
r=6
S=QQ[(symbol x)_0..(symbol x)_r];
L0=flatten apply(drop(gens S,1),g->({g-x_0,g+x_0})); L1=apply(toList(1..r-1),i->(x_i-x_(i+1))); L=join(L0,L1,{x_r-t*x_1,x_0}); A=arrangement L;
--The arrangement is free with exponents 1,3,3,3,3...
prune image der(A)

--create the restriction--
R=QQ[(symbol x)_1..(symbol x)_r];
L0=gens R; L1=apply(toList(1..r-1),i->(x_i-x_(i+1))); L=join(L0,L1,{x_r+x_1}); A=arrangement L;
--create multiplicity {2,2,..,2,1,...,1}
m=join(apply(r,i->2),apply(r,i->1))
--D(A) has pdim r-2
pdim image der A
--D(A,m) is free
pdim image der(A,m)

--------------------------------
--Example 5.7: Points in P^1----
--------------------------------

--first define the underlying polynomial ring
S=QQ[x,y];

--Define the B2 arrangement
t=-1;
L={x,y,x-y,x-t*y};
A=arrangement L

--Set a multiplicity
m={1,1,3,3}
D=derivationComplex(A,m);

--The 0th cohomology of D is the module of derivations
(prune HH^0 D)==(prune image der(A,m))

--change the parameter t from -1 to 2 to see
--the change in exponents of the module of derivations
prune HH^0 D

--Check that the formality complex S has the right form
(C,H)=formalityComplex(A)
--C is the formality complex (H encodes info about what flats
--the rows in the differentials of C correspond to)
C
--this is the coefficient matrix of A
C.dd_(0)
--this is the next differential
C.dd_(-1)

--Check that the complex J has the right form
J=kernelDerivationComplex(A,m)
J_(-2)

----------------------------------------
----Example 5.8: X3 arrangement --------
----------------------------------------

S=QQ[x,y,z];

--Define X3 arrangement--
t=-1;
L={x,y,z,x-t*y,x+z,y+z};
A=arrangement L;

--A non-free multiplicity--
m={2,2,2,2,2,2};
prune image der(A,m)

--can also check this isn't free by derivation complex
D=derivationComplex(A,m);
--First cohomology does not vanish, so not free
prune HH^1(D)

--A free multiplicity--
m={2,2,2,1,1,1};
--Free with exponents (3,3,3)--
prune image der(A,m)

--can also check freeness by derivation complex
D=derivationComplex(A,m);
--first cohomology vanishes, so free--
prune HH^1(D)

--Check formality complex
(C,H)=formalityComplex(A);
C.dd_0
C.dd_(-1)

--Check J complex--
J=kernelDerivationComplex(A,m)
--H^1 of derivation complex is H^2 of J complex
prune HH^2 J

------------------------------------------
---Example 6.5: Ziegler's Pair------------
------------------------------------------
S=QQ[x,y,z];

--six triple points of A1 lie on a smooth conic while six triple points of A2 do not
A1=arrangement({x,y,z,x+y+z,2*x+y+z,2*x+3*y+z,2*x+3*y+4*z,x+ 3*z,x+ 2*y+ 3*z});
A2=arrangement({x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z, 3*x+5*z, 3*x+4*y+5*z});

--A1 and A2 have the same intersection lattice--
L1=apply(flats A1,L->apply(L,f->f.flat));
L2=apply(flats A2,L->apply(L,f->f.flat));
L1==L2

--A1 is not formal but A2 is--
--Formality complexes change between A1 and A2 
--(there is a drop in rank from A2 to A1)
(C1,H1)=formalityComplex(A1);
C1
(C2,H2)=formalityComplex(A2);
C2

--A1 is not formal since HH^1 C2 does not vanish
prune HH^1 C1
--A2 is formal since HH^1 C2 vanishes
prune HH^1 C2

--non-freeness can be checked for small multiplicities (at most k)
--this will not finish if k is too big
k=2;
--create list of all multiplicities on A1 with largest multiplicity k
mList = apply(k,i->{i+1}); for i from 1 to 8 do(mList = flatten apply(k, j->(apply(mList,m->(append(m,j+1))))))
--loop through multiplicities and check that (A1,m) is not free
--On my machine, this takes about 60 seconds for k=2
scan(mList,m->(p = rank source presentation (prune image der(A1,m)); if p==0 then print m))