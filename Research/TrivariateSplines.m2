--Computations for the paper ``A lower bound for trivariate splines in large degree''
--loadPackage AlgebraicSplines --the function splineModule from the AlgebraicSplines package is used below

load("~/Dropbox/MN-Bounds/Macaulay2Scripts/BoundFunctions.m2") --(this function gives the local path to BoundFunctions.m2.  Alternatively,
--copy and paste the code from BoundFunctions.m2 into your Macaulay2 session)

--The commands below can be executed line by line in Macaulay2 to verify the computations in the paper 'A lower bound for trivariate splines in large degree'

----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Intro example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------

---Combinatorial Data to give the function 
F={{0,1,2,3},{0,1,2,7},{0,1,3,6},{0,2,3,5},{1,2,3,4},{0,1,6,7},{0,2,5,7},{0,3,5,6},{1,2,4,7},{1,3,4,6},{2,3,4,5},{0,5,6,7},{1,4,6,7},{2,4,5,7},{3,4,5,6}};

--define underlying polynomial ring (choose finite field for faster computation)
S=(ZZ/(nextPrime(1000)))[x,y,z,w]
--Geometric data for splineModule (this is symmetric)
V={{0,0,3},{1,1,-1},{-2,1,-1},{1,-2,-1},{0,0,-30},{-10,-10,10},{20,-10,10},{-10,20,10}};
--Move vertices to make generic
V'={{0,0,3},{1,1,-1},{-2,1,-1},{1,-2,-1},{-1,1,-30},{-11,-10,10},{20,-11,10},{-10,20,11}};

--Generic r=1 spline module
C1=splineModule(V',F,1,BaseRing=>S);
--Generic r=2 spline module
C2=splineModule(V',F,2,BaseRing=>S);
--Generic r=3 spline module (this takes a little while)
C3=splineModule(V',F,3,BaseRing=>S);
--Generic r=4 spline module (this takes a long time to compute!)
C4=splineModule(V',F,4,BaseRing=>S);

--r=1
--trivariateLower({1,d},F) computes our lower bound for this configuration when r=1 in degree d
--hilbertFunction(d,C1) computes the actual dimension for this configuration when r=1 in degree d
scan(toList(2..10),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({1,d},F),hilbertFunction(d,C1)}
	))
--r=2
--Dimension tests: the following code produces the data in Table 2
scan(toList(3..14),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({2,d},F),hilbertFunction(d,C2)}
	))
--r=3
--The following two blocks of code produces the data in Table 3
scan(toList(4..18),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({3,d},F),hilbertFunction(d,C3)}
	))
--r=4 (this takes a long time to compute!)
scan(toList(5..22),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({4,d},F),hilbertFunction(d,C4)}
	))

--Polynomial tests
--r=1
--trivariateLower({1},F) computes the polynomial form of our lower bound when r=1
--hilbertPolynomial(C1,Projective=>false) computes the actual polynomial giving the dimension of the spline
------space for large degree (the degree is recorded using the variable i)
trivariateLower({1},F)
hilbertPolynomial(C1,Projective=>false)
--r=2
--The polynomial for r=2 appears at the end of Section 2.1
trivariateLower({2},F)
hilbertPolynomial(C2,Projective=>false)
--The polynomials for r=3,4 appear in Section 5.1
--r=3
trivariateLower({3},F)
hilbertPolynomial(C3,Projective=>false)
--r=4 (this takes a long time to compute!)
trivariateLower({4},F)
hilbertPolynomial(C4,Projective=>false)

----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Section 5.3 Example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------

--Triangulation of 3D Torus
F={{0,3,7,15},{0,4,7,15},{0,4,12,15},{0,8,12,15},{0,8,11,15},{0,3,11,15},{0,1,4,12},{1,4,5,12},{1,5,12,13},{1,9,12,13},{1,8,9,12},{0,1,8,12},{1,2,5,13},{2,5,6,13},{2,6,13,14},{2,10,13,14},{2,9,10,13},{1,2,9,13},{2,3,6,14},{3,6,7,14},{3,7,14,15},{3,11,14,15},{3,10,11,14},{2,3,10,14}}

--define underlying polynomial ring (choose finite field for faster computation)
S=(ZZ/(nextPrime(1000)))[x,y,z,w]
--Geometric data for splineModule (this is symmetric)
V={{-1,1,0},{-1,-1,0},{1,-1,0},{1,1,0},{-2,2,0},{-2,-2,0},{2,-2,0},{2,2,0},{-1,1,1},{-1,-1,1},{1,-1,1},{1,1,1},{-2,2,1},{-2,-2,1},{2,-2,1},{2,2,1}};
--Move vertices to make generic
V'=apply(V,v->v+apply(3,i->random(1,10)/100));

--Generic r=1 spline module
C1=splineModule(V',F,1,BaseRing=>S);
--Generic r=2 spline module
C2=splineModule(V',F,2,BaseRing=>S);
--Generic r=3 spline module
C3=splineModule(V',F,3,BaseRing=>S);
--Generic r=4 spline module (this takes a long time to compute!)
C4=splineModule(V',F,4,BaseRing=>S);

--Dimension tests: the following produce that data in Table ??
--r=1
scan(toList(2..10),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({1,d},F),hilbertFunction(d,C1)}
	))
--r=2
scan(toList(3..14),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({2,d},F),hilbertFunction(d,C2)}
	))
--r=3
scan(toList(4..18),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({3,d},F),hilbertFunction(d,C3)}
	))
--r=4 (this takes a long time to compute!)
scan(toList(5..22),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({4,d},F),hilbertFunction(d,C4)}
	))

--Polynomial tests
--r=1
trivariateLower({1},F)
hilbertPolynomial(C1,Projective=>false)
--r=2
trivariateLower({2},F)
hilbertPolynomial(C2,Projective=>false)
--r=3
trivariateLower({3},F)
hilbertPolynomial(C3,Projective=>false)
--r=4 (this takes a long time to compute!)
trivariateLower({4},F)
hilbertPolynomial(C4,Projective=>false)


----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Section 5.2 Example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------

--3D Morgan Scot with cavity
--Combinatorial data
F={{0,1,2,7},{0,1,3,6},{0,2,3,5},{1,2,3,4},{0,1,6,7},{0,2,5,7},{0,3,5,6},{1,2,4,7},{1,3,4,6},{2,3,4,5},{0,5,6,7},{1,4,6,7},{2,4,5,7},{3,4,5,6}};

--define underlying polynomial ring (choose finite field for faster computation)
S=(ZZ/(nextPrime(1000)))[x,y,z,w]
--Geometric data for splineModule (this is symmetric)
V={{0,0,3},{1,1,-1},{-2,1,-1},{1,-2,-1},{0,0,-30},{-10,-10,10},{20,-10,10},{-10,20,10}};
--Move vertices to make generic
V'={{0,0,3},{1,1,-1},{-2,1,-1},{1,-2,-1},{-1,1,-30},{-11,-10,10},{20,-11,10},{-10,20,11}};

--Generic r=1 spline module
C1=splineModule(V',F,1,BaseRing=>S);
--Generic r=2 spline module
C2=splineModule(V',F,2,BaseRing=>S);
--Generic r=3 spline module
C3=splineModule(V',F,3,BaseRing=>S);
--Generic r=4 spline module (this takes a long time to compute!)
C4=splineModule(V',F,4,BaseRing=>S);

--Dimension tests: the following produce that data in Table ??
--r=1
scan(toList(2..10),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({1,d},F),hilbertFunction(d,C1)}
	))
--r=2
scan(toList(3..14),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({2,d},F),hilbertFunction(d,C2)}
	))
--r=3
scan(toList(4..18),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({3,d},F),hilbertFunction(d,C3)}
	))
--r=4 (this takes a long time to compute!)
scan(toList(5..22),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({4,d},F),hilbertFunction(d,C4)}
	))

--Polynomial tests
--r=1
trivariateLower({1},F)
hilbertPolynomial(C1,Projective=>false)
--r=2
trivariateLower({2},F)
hilbertPolynomial(C2,Projective=>false)
--r=3
trivariateLower({3},F)
hilbertPolynomial(C3,Projective=>false)
--r=4 (this takes a long time to compute!)
trivariateLower({4},F)
hilbertPolynomial(C4,Projective=>false)

----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Section 5.4 Example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------

---Cube Octahedron (Non-simplicial)
--Combinatorial data (need to give it explicitly since its a non-simplicial example)
F = {{0, 1, 2, 3, 4, 5}, {0, 8, 9, 12, 13}, {1, 6, 7, 10, 11}, {2, 7, 8, 11, 12}, {3, 6, 9, 10, 13}, {4, 10, 11, 12, 13}, {5, 6, 7, 8, 9}, {0, 2, 8, 12}, {0, 3, 9, 13}, {0, 4, 12, 13}, {0, 5, 8, 9}, {1, 2, 7, 11}, {1, 3, 6, 10}, {1, 4, 10, 11}, {1, 5, 6, 7}, {2, 4, 11, 12}, {3, 4, 10, 13}, {3, 5, 6, 9}, {2, 5, 7, 8}, {0, 2, 4, 12}, {0, 2, 5, 8}, {0, 3, 4, 13}, {0, 3, 5, 9}, {1, 2, 4, 11}, {1, 2, 5, 7}, {1, 3, 4, 10}, {1, 3, 5, 6}};
numFaces = 56 --number of 2-faces
eulerC = 1 --Euler characteristic is 1
edgeIncidences = apply(36,i->4) --each of the 36 interior edges is surrounded by 4 two-faces
intStarIncidences = apply(6,i->{16,{4,4,4,4,4,4,4,4},9}) --each interior vertex is contained in 16 two-faces, 8 edges which are all surrounded by 4 two-faces, and 16 polytopes
bndryStarIncidences = apply(8,i->{9,{4,4,4},7}) --each boundary vertex is contained in 9 interior two-faces, 3 interior edges each of which is surrounded by 4 two-faces, and 7 polytopes

--Geometric data
--Define polynomial ring over smaller finite field to make computations reasonable
S=(ZZ/nextPrime(1000))[x,y,z,w]
--Vertex coordinates (symmetric case)
V = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {-2, -2, -2}, {-2, 2, -2}, {2, 2, -2}, {2, -2, -2}, {-2, -2, 2}, {-2, 2, 2}, {2, 2, 2}, {2, -2, 2}};
--Wiggle vertices to make this generic
V'=apply(V,v->apply(v,i->i+random(1,10)/1000))

--These spline module computations take quite some time, and possibly won't finish for r=3,4
--Generic r=1 spline module (this takes a long time)
C1=splineModule(V',F,1,BaseRing=>S);
--Generic r=2 spline module
C2=splineModule(V',F,2,BaseRing=>S);
--Generic r=3 spline module (this computation takes a long time!)
C3=splineModule(V',F,3,BaseRing=>S);
--Generic r=4 spline module (this computation may not finish)
C4=splineModule(V',F,4,BaseRing=>S);

--Dimension tests: the following produce the data in Table 
--r=1
scan(toList(2..10),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({1,numFaces,eulerC,d},edgeIncidences,intStarIncidences,bndryStarIncidences),hilbertFunction(d,C1)}
	))
--r=2
scan(toList(3..14),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({2,numFaces,eulerC,d},edgeIncidences,intStarIncidences,bndryStarIncidences),hilbertFunction(d,C2)}
	))
--r=3
scan(toList(4..18),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({3,numFaces,eulerC,d},edgeIncidences,intStarIncidences,bndryStarIncidences),hilbertFunction(d,C3)}
	))
--r=4
scan(toList(5..22),d->(
	<<endl;
	print{d,binom(d+3,3),trivariateLower({4,numFaces,eulerC,d},edgeIncidences,intStarIncidences,bndryStarIncidences),hilbertFunction(d,C4)}
	))

--Polynomial tests
--r=1
trivariateLower({1,numFaces,eulerC},edgeIncidences,intStarIncidences,bndryStarIncidences)
hilbertPolynomial(C1,Projective=>false)
--r=2
trivariateLower({2,numFaces,eulerC},edgeIncidences,intStarIncidences,bndryStarIncidences)
hilbertPolynomial(C2,Projective=>false)
--r=3
trivariateLower({3,numFaces,eulerC},edgeIncidences,intStarIncidences,bndryStarIncidences)
hilbertPolynomial(C3,Projective=>false)
--r=4
trivariateLower({4,numFaces,eulerC},edgeIncidences,intStarIncidences,bndryStarIncidences)
hilbertPolynomial(C4,Projective=>false)

LBF:=d->(7*binom(d+3,3)-52*binom(d+1,3)+72*binom(d,3)-12)

scan(25,d->(
	<<endl;
	print{d,binom(d+3,3),LBF(d),genericPredictedLower(r,d,faces,edgeIncidences,ivertexIncidences,evertexIncidences),hilbertFunction(d,C)}
	))

scan(25,d->(
	<<endl;
	print{d,binom(d+3,3),LBF(d)-hilbertFunction(d,C)}
	))
