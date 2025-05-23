--this file can be loaded into a Macaulay2 session and will give access to the functions it defines

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
    
    
isThreeGenerated = method()
isThreeGenerated = A -> (
    R := ring A;
    K := coefficientRing(R);
    M := coefficients A;
    RELS = ker M;
    ncM := numcols M;
    nrM := numrows M;
    FL2 := apply(flats(2,A),f->toList f);--get codimension two flats, turn into lists
    FL2trip := select(FL2,f->#f>=3);
    --If there are no flats of codimension two, then there are no relations of length three.
    if #FL2trip == 0 then(
	if (rank RELS)>0 then(
	    return false
	    )else(
	    return true
	    ));
    --The next function gets the length three relations as a submodule of all relations
    LCLRELS := sum apply(FL2trip, f->(
	    fsz := #f;
	    idsz := ncM-fsz;
	    c := 0;
	    i := 0;
	    ID := id_(K^(ncM-fsz));
	    relMat := {};
	    while i+c<ncM do(
		if isMember(i+c,f) then(
		    relMat = append(relMat,join(flatten entries(M_(i+c)),apply(idsz,i->0)));
		    c = c+1
		    )else(
		    relMat = append(relMat,join(apply(nrM,i->0),flatten entries(ID_i)));
		    i = i+1
		    )
		);
	    ker transpose matrix relMat
	    ));
    QUOT := prune(RELS/LCLRELS);
    if (rank QUOT)>0 then return false else true
    )
    
	    
