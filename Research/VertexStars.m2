--Computations for the paper ``A lower bound for splines on tetrahedral vertex stars''
loadPackage AlgebraicSplines --the function splineModule from the AlgebraicSplines package is used below

load("~/Macaulay2Scripts/BoundFunctions.m2") --(this function gives the local path to BoundFunctions.m2.  Alternatively,
--copy and paste the code from BoundFunctions.m2 into your Macaulay2 session)

----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Section 6.1 example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------

--Generic bipyramid over a pentagon
--Triangles listed are the triangles of the link (this is the input for the function
--"closedVertexStarLower"
F={{1,2,6},{2,3,6},{3,4,6},{4,5,6},{1,5,6},{1,2,7},{2,3,7},{3,4,7},{4,5,7},{1,5,7}};

--define underlying polynomial ring
S=QQ[x,y,z];
--vertex coordinates for non-central vertices, used to compute actual dimension
V={{1,0,0},{1,4,0},{-1,1,0},{-2,-1,0},{1,-2,0},{0,0,1},{-2,1,-2}};
--make the coordinates generic by adding a small random vector to each vertex
V'=prepend({0,0,0},V+apply(V,v->v+apply(3,i->random(1,10)/100)));
--For input into the function "splineModule", we need to append the central vertex to each triangle
F1={{0,1,2,6},{0,2,3,6},{0,3,4,6},{0,4,5,6},{0,1,5,6},{0,1,2,7},{0,2,3,7},{0,3,4,7},{0,4,5,7},{0,1,5,7}}
--r=1 spline module
C1=splineModule(V',F1,1,Homogenize=>false,BaseRing=>S);
--r=2 spline module
C2=splineModule(V',F1,2,Homogenize=>false,BaseRing=>S);
--r=3 spline module (takes a few seconds)
C3=splineModule(V',F1,3,Homogenize=>false,BaseRing=>S);
--r=4 spline module (takes 10 seconds or so)
C4=splineModule(V',F1,4,Homogenize=>false,BaseRing=>S);

--the following produces the output for the left portion of Table 1
--r=1
scan(toList(2..9),d->print {d,binom(d+2,2), closedVertexStarLower(1,d,F), hilbertFunction(d,C1)})
--r=2
scan(toList(3..11),d->print {d,binom(d+2,2), closedVertexStarLower(2,d,F), hilbertFunction(d,C2)})
--r=3
scan(toList(4..12),d->print {d,binom(d+2,2), closedVertexStarLower(3,d,F), hilbertFunction(d,C3)})
--r=4
scan(toList(5..13),d->print {d,binom(d+2,2), closedVertexStarLower(4,d,F), hilbertFunction(d,C4)})


----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Section 6.2 example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
--Non-generic bipyramid over a pentagon

--combinatorial data for non-generic bipyramid
numfaces=15 --15 interior faces
numslopeslist={5,5,3,3,3,3,3} --the five rays on z=0 have only 3 slopes

--Tetrahedra of bipyramid
F1={{0,1,2,6},{0,2,3,6},{0,3,4,6},{0,4,5,6},{0,1,5,6},{0,1,2,7},{0,2,3,7},{0,3,4,7},{0,4,5,7},{0,1,5,7}}
--Need coordinates for exact computation
V={{1,0,0},{1,4,0},{-1,1,0},{-2,-1,0},{1,-2,0},{0,0,1},{-2,1,-2}};
--apply random small translations parallel to xy plane
V'=prepend({0,0,0},V+apply(V,v->v+append(apply(2,i->random(1,10)/100),0)));
--r=1 spline module
C1ng=splineModule(V',F1,1,Homogenize=>false,BaseRing=>S);
--r=2 spline module
C2ng=splineModule(V',F1,2,Homogenize=>false,BaseRing=>S);
--r=3 spline module
C3ng=splineModule(V',F1,3,Homogenize=>false,BaseRing=>S);
--r=4 spline module
C4ng=splineModule(V',F1,4,Homogenize=>false,BaseRing=>S);

--the following produces the output for the right portion of Table 1
--r=1
scan(toList(2..9),d->print {d,binom(d+2,2)+binom(d,2), closedVertexStarLower(1,d,numfaces,numslopeslist), hilbertFunction(d,C1ng)})
--r=2
scan(toList(3..11),d->print {d,binom(d+2,2)+binom(d-1,2), closedVertexStarLower(2,d,numfaces,numslopeslist), hilbertFunction(d,C2ng)})
--r=3
scan(toList(4..12),d->print {d,binom(d+2,2)+binom(d-2,2), closedVertexStarLower(3,d,numfaces,numslopeslist), hilbertFunction(d,C3ng)})
--r=4
scan(toList(5..13),d->print {d,binom(d+2,2)+binom(d-3,2), closedVertexStarLower(4,d,numfaces,numslopeslist), hilbertFunction(d,C4ng)})

----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Section 6.3 example
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------

--combinatorial data for generic cube
numfaces=12 --12 interior faces
numslopeslist=apply(8,i->3) --all edges have three incident slopes

--geometric data for generic cube
--coordinates for symmetric cube
V={{1,-1,1},{-1,-1,1},{-1,1,1},{1,1,1},{1,-1,-1},{-1,-1,-1},{-1,1,-1},{1,1,-1}};
--add small vector to each of vertices to get generic coordinates
V'=prepend({0,0,0},apply(V,L->apply(L,i->i+random(-1,1)/10)))
--Faces for generic cube
F={{0,1,2,3,4},{0,1,2,5,6},{0,1,4,5,8},{0,2,3,6,7},{0,3,4,7,8},{0,5,6,7,8}};

S=QQ[x,y,z]
--Spline modules
--r=1
S1=splineModule(V',F,1,Homogenize=>false,BaseRing=>S);
--r=2
S2=splineModule(V',F,2,Homogenize=>false,BaseRing=>S); --takes a while
--r=3
S3=splineModule(V',F,3,Homogenize=>false,BaseRing=>S); --takes a very long time, replace QQ by finite field for faster computation

--the following produces the output for Table 2
--r=1
scan(toList(2..9),d->print {d,binom(d+2,2), closedVertexStarLower(1,d,numfaces,numslopeslist), hilbertFunction(d,S1)})
--r=2
scan(toList(3..11),d->print {d,binom(d+2,2), closedVertexStarLower(2,d,numfaces,numslopeslist), hilbertFunction(d,S2)})
--r=3
scan(toList(5..15),d->print {d,binom(d+2,2), closedVertexStarLower(3,d,numfaces,numslopeslist), hilbertFunction(d,S3)})

----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
---------Additional examples (forthcoming)
----------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------
