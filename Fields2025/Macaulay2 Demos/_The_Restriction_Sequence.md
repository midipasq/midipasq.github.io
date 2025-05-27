# A deletion-restriction sequence

## Loading functions and packages

Upload "Scripts_for_splines_on_hyperplanes.m2" or load the relevant packages/functions below.

```
loadPackage("AlgebraicSplines",Reload=>true)
loadPackage("HyperplaneArrangements",Reload=>true)

zonotopeEdgeLabels = method()
--Inputs: A: a hyperplane arrangement A
--Outputs: A sequence (E,I) where
--        E records the edges between vertices of Z(A)
--        I records the edge labels using the smoothness distribution r
zonotopeEdgeLabels = A ->(
    S := ring A;--get the ring of A
    hyps := A#(first keys A); --get the list of hyperplanes definining A
    M := coefficients A; --get the coefficient matrix of A
    normals := entries transpose M;
    Z := sum apply(numcols M,j->convexHull((-M_{j})|M_{j}));--get the zonotope Z(A)
    V := entries transpose vertices Z; --get vertices of Z(A)
    --the remaining code inductively builds the list of edges and edge labels
    --by cycling through the vertices of Z(A) as recorded as columns of the matrix "vertices Z"
    --and adding the vector 2v for each vector in 'normals'.  If the result is a vertex of Z(A),
    --then the corresponding edge (recorded as a pair {i,j} of column indices of the matrix "vertices Z")
    --is added to E and the corresponding edge label is added to I.
    E := {};
    I := {};
    vtxcnt := 0;
    while vtxcnt<(#V) do(
	nrmcnt := 0;
	vtx := V_(vtxcnt);
	while nrmcnt<#normals do(
	    nrm := 2*normals_nrmcnt;
	    vtxc := vtx+nrm;
	    p := position(V,v->v==vtxc);
	    if not p===null then(
		E = append(E,{vtxcnt,p});
		I = append(I,hyps_(nrmcnt));
		);
	    nrmcnt=nrmcnt+1;
	    );
	vtxcnt = vtxcnt+1
	);
    (E,I,Z)
    )

assignSmoothness = method()
assignSmoothness = (A,r,I)->(
    hyps := A#(first keys A);
    smthHash := hashTable apply(length r,i->{hyps_i,r_i});
    apply(I,h->h^(smthHash#h+1))
    )

```

## Deletion

```
S=QQ[x,y]
A=arrangement{x,y,x-y,x+y}
(E,L,Z)=zonotopeEdgeLabels A
r={2,2,3,3}
I=assignSmoothness(A,r,L)
GS=generalizedSplines(E,I)
prune GS
```
```
rd={2,2,3,2}
Id=assignSmoothness(A,rd,L)
GSd=generalizedSplines(E,Id)
prune GSd
