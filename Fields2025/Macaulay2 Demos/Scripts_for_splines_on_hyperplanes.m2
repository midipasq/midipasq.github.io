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
    L := {};
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
		L = append(L,hyps_(nrmcnt));
		);
	    nrmcnt=nrmcnt+1;
	    );
	vtxcnt = vtxcnt+1
	);
    (E,L,Z)
    )

assignSmoothness = method()
--Inputs:   A: hyperplane arrangement
--	    r: smoothness distribution (same order as hyperplanes)
--	    L: linear forms labeling edges
--Outputs:  I: linear forms of L raised to powers assigned by r
assignSmoothness = (A,r,L)->(
    hyps := A#(first keys A);
    smthHash := hashTable apply(length r,i->{hyps_i,r_i});
    apply(L,h->h^(smthHash#h+1))
    )
    
    
isThreeGenerated = method()
--Inputs: A: Hyperplane arrangement
--Outputs: Boolean telling whether A is three-generated or not
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
    
sageLabels = method()
sageLabels = (Z,E,I)->(
    V := entries transpose vertices Z;
    V1 := toString V;
    V2 := replace("[{]","[",V1);
    V3 := replace("[}]","]",V2);
    V4 := "V="|V3;
    E = apply(E,e->sort(e));
    Es0 := toString E;
    Es1 := Es0_(1,length(Es0)-2);
    Es2 := replace("[{]","(",Es1);
    Es3 := replace("[}]",")",Es2);
    Es4 := "E="|"["|Es3|"]";
    Is := concatenate apply(I,f->"'"|toString(f)|"',");
    Is1 := "["|(Is_(0,length(Is)-1))|"]";
    Is2 := "edgelabels="|Is1;
    print V4;
    print Es4;
    print Is2
    )
    
localSyzygyModule = method()
localSyzygyModule = (E,I,Z) -> (
    S := ring first I;
    r := flatten apply(I,f->degree f);
    nhyps := #I;
    FL20 := apply(faces(1,Z),f->first f);--get codim 1 faces, turn into list
    FL2 := apply(FL20,f->positions(E,e->isSubset(e,f))); --turn into list of edges
    --The next function gets the local syzygies as a submodule of the appropriate free module
    sum apply(FL2, f->(
	    fsz := #f;
	    idsz := nhyps-fsz;
	    c := 0;
	    i := 0;
	    ID := id_(S^(nhyps-fsz));
	    relMat := {};
	    while i+c<nhyps do(
		if isMember(i+c,f) then(
		    relMat = append(relMat,join({I_(i+c)},apply(idsz,i->0_S)));
		    c = c+1
		    )else(
		    relMat = append(relMat,join({0_S},flatten entries(ID_i)));
		    i = i+1
		    )
		);
	    phi=map(S^(1+idsz),S^(-r),transpose matrix relMat);
	    ker phi
	    ))
    )

directSumLocalSyzygyModule = (E,I,Z) -> (
    S := ring first I;
    FL20 := apply(faces(1,Z),f->first f);--get codim 1 faces, turn into list
    FL2 := apply(FL20,f->positions(E,e->isSubset(e,f))); --turn into list of edges
    --The next function gets the local syzygies as a submodule of the appropriate free module
    image directSum apply(FL2, f->(
	    syz matrix {I_f}
	    ))
    )


reducedLocalSyzygyModule = method()
reducedLocalSyzygyModule = (A,r)->(
    S := ring first L;
    hyps := A#(first keys A);
    nhyps := #hyps;
    FL2 := apply(flats(2,A),f->toList f);--get codimension two flats, turn into lists
    --The next function gets the local syzygies as a submodule of the appropriate free module
    sum apply(FL2, f->(
	    fsz := #f;
	    idsz := nhyps-fsz;
	    c := 0;
	    i := 0;
	    ID := id_(S^(nhyps-fsz));
	    relMat := {};
	    while i+c<nhyps do(
		if isMember(i+c,f) then(
		    relMat = append(relMat,join({hyps_(i+c)^(r_(i+c)+1)},apply(idsz,i->0_S)));
		    c = c+1
		    )else(
		    relMat = append(relMat,join({0_S},flatten entries(ID_i)));
		    i = i+1
		    )
		);
	    phi=map(S^(1+idsz),S^(apply(r,i->-i-1)),transpose matrix relMat);
	    ker phi
	    ))
    )
