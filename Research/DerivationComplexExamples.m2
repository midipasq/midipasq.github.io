--File for testing examples of chain complexes in the file DerivationComplex.m2--
--package "HyperplaneArrangements" is automatically loaded in the following file
load "DerivationComplex.m2" --this line should be modified to the path where the file is located, for instance load "~/Documents/M2/DerivationComplex.m2"

--Let A be an arrangement with n hyperplanes in a vector space V of dimension l.
--Set S=sym(V*).  The function formalityComplex(A) returns a sequence (C,H). Ignore H
--for now.  C is a chain complex which satisfies the following property:
--If A is k-formal but not (k+1)-formal, then HH^i(C)=0 for i=1,..,k-1
--and HH^k(C) is not zero.  The chain complex is the dual of a chain complex
--in a paper of Stefan Tohaneanu; the definition of k-formal appears in a paper
--of Brandt and Terao.
S=QQ[w,x,y,z]
A=typeA(3,S)
(C,H)=formalityComplex(A);
--The A3 arrangement is k-formal for all k, but not essential, so HH^0 is not zero but all higher cohomologies vanish.
scan(length C+1,i->print({i,prune HH^i C}))
--taking HH^i is same as taking HH_(-i) (this is why C is indexed
--negatively
prune HH C

--Here is an example of Brandt and Terao which is formal but not three-formal.
S=QQ[x,y,z,w];
L={z,z-w,y,y+z-2*w,x,x+z-2*w,y+2*z-2*w,x+2*z-2*w,x+y+z-2*w,w};
A=arrangement L;
--this takes some time because it calls the flats method
(C,H)=formalityComplex(A);
--A is essential and formal (i.e. 2-formal), but not 3-formal, so H^0=H^1=0 but H^2 does not vanish
prune HH^0 C
prune HH^1 C
prune HH^2 C

--Here is an example of Yuzvinsky/Ziegler which shows formality is not combinatorial--
S=QQ[x,y,z];
--six triple points of A1 lie on a smooth conic while six triple points of A2 do not
A1=arrangement({x,y,z,x+y+z,2*x+y+z,2*x+3*y+z,2*x+3*y+4*z,x+ 3*z,x+ 2*y+ 3*z});
A2=arrangement({x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z, 3*x+5*z, 3*x+4*y+5*z});
--A1 and A2 have the same intersection lattice--
L1=apply(flats A1,L->apply(L,f->f.flat))
L2=apply(flats A2,L->apply(L,f->f.flat))
L1==L2
--A1 is not formal but A2 is--
(C1,H1)=formalityComplex(A1);
(C2,H2)=formalityComplex(A2);
--A1 is not formal since HH^1 C2 does not vanish
prune HH^1 C1
--A2 is formal since HH^1 C2 vanishes
prune HH^1 C2

--Let m be a list of positive integers, multiplicities on the hyperplanes of A
--The function derivationComplex(A,m) returns a chain complex D with the following
--properties:
--1. D is a quotient of the chain complex returned by formalityComplex(A)
--2. HH^0(D) is the module of multi-derivations D(A,m)
--3. HH^i(D)=0 for i>0 if and only if D(A,m) is free
S=QQ[w,x,y,z];
A=typeA(3,S);
D=derivationComplex(A)
--A is free with exponents 0,1,2,3
prune HH^0 D
--higher cohomologies vanish
prune HH D
--a non-free multiplicity on A3
m={1,1,2,2,1,1};
D=derivationComplex(A,m);
prune HH^0 D
--We can also tell (A,m) is not free because HH^1 D does not vanish
prune HH^1 D

--the following function produces the chain complex which is the kernel
--of the surjection formalityComplex(A)->derivationComplex(A,m)
S=QQ[w,x,y,z];
A=typeA(3,S);
kD=kernelDerivationComplex(A);
(C,H)=formalityComplex(A);
D=derivationComplex(A);
D==coker(inducedMap(C,kD))
--HH^1 kD is the derivations D(A,m) without those of degree zero
prune HH^1 kD
prune HH^0 D

--Here are the details of the sequence (C,H) which appears as output of formalityComplex(A).
--C is a chain complex, C_0=S^l, C_(-1)=S^n, and the rows of C.dd_(0) are the coefficients
--of forms defining hyperplanes of A, up to a global sign change (i.e. the transpose of coefficients(A)).
S=QQ[w,x,y,z]
A=typeA(3,S)
(C,H)=formalityComplex(A);
C_0
C_(-1)
C.dd_0
--C_(-2) has a basis in 1-1 correspondence with relations among three forms of A (i.e. these correspond
--to relations around codimension 2 flats where at least 3 hyperplanes meet), and the rows of C.dd_(-1)
--are these relations.
C.dd_(-1)
--Iteratively, the rows of C_(-d) are relations on rows of C_(-d+1) which correspond
--to flats of codimension d where at least d+1 hyperplanes meet.
--H is a hash table.  H_d encodes the flats where non-trivial relations occur on rows
--of C.dd_(-d), and which rows of C.dd_(-d) correspond to this flat.  This information
--is recorded as Flat=>List, where List records the indices of rows of C.dd_(-d) that
--have a non-trivial relation around the given Flat.
H_0
H_1
H_2