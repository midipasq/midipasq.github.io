--Functions to compute bounds in papers:
---"A lower bound for splines on tetrahedral vertex stars"
---"A lower bound on trivariate splines in large degree"

-----------------------------------------------------------------------
--redefine the binomial function to fit Hilbert function convention
-----------------------------------------------------------------------
binom=(n,k)->(
    if n<k then 0 else binomial(n,k)
    )
-----------------------------------------------------------------------

edgeTerm = method(Options=>{
	symbol BaseRing => null
	}
    )
-----------------------------------------------------------------------
--The main edge term (dimension of syz J(tau)_d)
--Inputs:
----r: degree of smoothness
----d: polynomial degree
----ngens: number of generators of J(tau) (minimum of r+2 and number of slopes
-------meeting at tau)
----ambdim: Either two (if computing for vertex stars) or three (if computing
-------for tetrahedral complexes)
-----------------------------------------------------------------------
edgeTerm(ZZ,ZZ,ZZ,ZZ) := opts-> (r,d,ngens,ambdim)->(
    if ngens<=1 then(
	0
	)else if ngens==2 then(
	binom(d+ambdim-2*(r+1),ambdim)
	)else(
	q := floor(ngens*(r+1)/(ngens-1));
	a := ngens*(r+1)-q*(ngens-1);
	b := ngens-1-a;
	a*binom(d+ambdim-q-1,ambdim)+b*binom(d+ambdim-q,ambdim)
	)
    )

edgeTerm(ZZ,ZZ,ZZ) := opts->(r,ngens,ambdim)->(
    if (class opts.BaseRing) === PolynomialRing then(
	R := opts.BaseRing
	)else(
	R = QQ[symbol d];
	);
    if ngens<=1 then(
	return 0_R
	)else if ngens==2 then(
	return binomial(d+ambdim-2*(r+1),ambdim)
	)else(
	q := floor(ngens*(r+1)/(ngens-1));
	a := ngens*(r+1)-q*(ngens-1);
	b := ngens-1-a;
	return a*binomial(d+ambdim-q+1,ambdim)+b*binomial(d+ambdim-q,ambdim)
	)
    )


openVertexStarLower=method()
-----------------------------------------------------------------------
--Lower bound for open vertex star (Schumaker's lower bound)
--Inputs:
---r: smoothness
---d: degree to compute lower bound in
---Additional data depending on input type for the open vertex star (see below)
-----------------------------------------------------------------------

openVertexStarLower(ZZ,ZZ,ZZ,List) := (r,d,numfaces,numslopeslist)->(
-----------------------------------------------------------------------
--Inputs:
---r: smoothness
---d: degree to compute lower bound in
---numfaces: number of interior dimension two faces in the open vertex star
---numslopeslist: list whose entries record number of 2-dim faces
----- meeting at each interior one-dim face (or the number of distinct planes
----- spanned by these 2-dim faces)
--Outputs:
---A lower bound for dimension of homogeneous splines on the open vertex star
----(exact in high degree if number of slopes is given for each vertex)
-----------------------------------------------------------------------
    --the global splines
    glbTerm := binom(d+2,2);
    --replace numslopes for each ray with min of numslopes and r+2
    ntlist := apply(numslopeslist,s->min(r+2,s));
    --the degree r+1 contribution
    faceTerm := (numfaces-sum(ntlist))*binom(d+1-r,2);
    --the contribution from each edge
    edgeT := sum apply(ntlist, s-> edgeTerm(r,d,s,2));
    glbTerm+faceTerm+edgeT
    )


openVertexStarLower(ZZ,ZZ,List) := (r,d,F)->(
-----------------------------------------------------------------------
--Inputs:
---r: smoothness
---d: degree to compute lower bound in
---F: list of faces of Delta according to the vertices they contain (omit
------- the central vertex from each face!)  Must be simplicial.
--Output: A lower bound for dimension of homogeneous splines on the open vertex star
----(exact in high degree only if vertex coordinates are generic)
-----------------------------------------------------------------------
    --Check to see that input is simplicial
    if all(F,s->#s==3)then(
	--get the rays (dim 1 faces)
	rays := unique flatten F;
	--get all dim 2 faces
	allFaces := unique apply(flatten apply(F,f->subsets(f,2)),e->sort(e));
	--get interior dim 2 faces
	intFaces := select(allFaces, e->#select(F,f->isSubset(set e,set f))==2);
	--get number of interior dim 2 faces
	numfaces := length intFaces;
	--get boundary dim 2 faces
	bfaces := select(allFaces, e->#select(F,f->isSubset(set e,set f))==1);
	--get boundary rays
	brays := unique flatten bfaces;
	--get interior rays
	intrays := rays-set(brays);
	--plug this input into the previously defined function
	numslopeslist := apply(intrays,r->#select(intFaces,e->member(r,e)));
	return openVertexStarLower(r,d,numfaces,numslopeslist)
	)else(
	return error "Must be simplicial for this form of input"
	)
    )

closedVertexStarLower=method()
-----------------------------------------------------------------------
--Lower bound for closed vertex star
--Inputs:
---r: smoothness
---d: degree to compute lower bound in
---Additional data depending on input type for the closed vertex star (see below)
-----------------------------------------------------------------------

closedVertexStarLower(ZZ,ZZ,ZZ,List) := (r,d,numfaces,numslopeslist)->(
-----------------------------------------------------------------------
--Inputs:
---r: smoothness
---d: degree to compute lower bound in
---numfaces: number of interior dimension two faces in the open vertex star
---numslopeslist: list whose entries record number of 2-dim faces
----- meeting at each interior one-dim face (or the number of distinct planes
----- spanned by these 2-dim faces)
--Outputs:
---A lower bound for dimension of homogeneous splines on the closed vertex star
----(exact in high degree if number of slopes is given for each vertex)
-----------------------------------------------------------------------
    --the global splines
    glbTerm := binom(d+2,2);
    --replace numslopes for each ray with min of numslopes and r+2
    ntlist := apply(numslopeslist,s->min(r+2,s));
    --the degree r+1 contribution
    faceTerm := (numfaces-sum(ntlist))*binom(d+1-r,2);
    --the contribution from each edge
    edgeT := sum apply(ntlist, s-> edgeTerm(r,d,s,2));
    2*glbTerm+faceTerm+edgeT
    )


closedVertexStarLower(ZZ,ZZ,List) := (r,d,F)->(
-----------------------------------------------------------------------
--Inputs:
---r: smoothness
---d: degree to compute lower bound in
---F: list of faces of Delta according to the vertices they contain (omit
------- the central vertex from each face!)  Must be simplicial.
--Output: A lower bound for dimension of homogeneous splines on the closed vertex star
----(exact in high degree only if vertex coordinates are generic)
-----------------------------------------------------------------------
    --Check to see that input is simplicial
    if all(F,s->#s==3)then(
	--get the rays (dim 1 faces)
	rays := unique flatten F;
	--get all dim 2 faces
	allFaces := unique apply(flatten apply(F,f->subsets(f,2)),e->sort(e));
	--get interior dim 2 faces
	numfaces := length allFaces;
	--plug this input into the previously defined function
	numslopeslist := apply(rays,r->#select(allFaces,e->member(r,e)));
	return closedVertexStarLower(r,d,numfaces,numslopeslist)
	)else(
	return error "Must be simplicial for this form of input"
	)
    )

vertexContribution=method()
-----------------------------------------------------------------------
--Compute contribution to trivariate lower bound from each vertex
--Inputs: 
---r,numfaces, and numslopeslist which are same as for openVertexStarLower
-----and closedVertexStarLower
---numpolytopes: integer, the number of polytopes in case the vertex 
-----star is open, zero if its closed
-----------------------------------------------------------------------
vertexContribution(ZZ,ZZ,List,ZZ) := (r,numfaces,numslopeslist,numpolytopes)->(
    --Initialize vertex contribution
    Ng := 0;
    if numpolytopes==0 then(
	--Define Dgamma
	if numfaces==4 then(
	    Dg := 2*r
	    );
	if numfaces==5 then(
	    Dg = floor((5*r+2)/3)
	    );
	if numfaces>5 then(
	    Dg = floor((3*r+1)/2)
	    );
	--Initialize vertex contribution for closed star
	Ng = -binom(r+3,3);
	j := r+1;
	--Allow additional negative contributions when j<= Dg
	while j<= Dg do(
	    Ng = Ng+binom(j+2,2)-closedVertexStarLower(r,j,numfaces,numslopeslist);
	    j = j+1
	    );
	--get number of interior edges
	numrays := length numslopeslist;
	--Maximum of binom(d+2,2)-closedVertexStarLower(r,d,numfaces,numslopeslist) occurs at apex
	apex := (r+1)*(numfaces/(1-numrays+numfaces)-3/2);
	flag := 0;
	--Add any additional contributions to Ng until they are no longer positive
	while flag == 0  do(
	    toAdd := binom(j+2,2)-closedVertexStarLower(r,j,numfaces,numslopeslist);
	    if toAdd>0 then(
		Ng = Ng+toAdd;
		j = j+1
		)else(
		if j>apex then(
		    flag = 1
		    )else(
		    j = j+1
		)))
    )else(
    --Maximum of binom(d+2,2)-closedVertexStarLower(r,d,numfaces,numslopeslist) occurs at apex
    apex = (r+1)*(numfaces/(numpolytopes-1)-3/2);
    flag = 0;
    j = r+1;
    --Add any additional contributions to Ng until they are no longer positive
    while flag == 0  do(
	toAdd = binom(j+2,2)-openVertexStarLower(r,j,numfaces,numslopeslist);
	if toAdd>0 then(
	    Ng = Ng+toAdd;
	    j = j+1
	    )else(
	    if j>apex then(
		flag = 1
		)else(
		j = j+1
		))));
    Ng
    )
    

trivariateLower=method()
-----------------------------------------------------------------------
--Compute lower bound in large degree
--Inputs: depends on input type for the tetrahedral complex (see below)
--Output: lower bound in large degree for splines on the tetrahedral complex
-----------------------------------------------------------------------

trivariateLower(List,List,List,List) := (rdfe,edgeIncidences,intStarIncidences,bndryStarIncidences)->(
-----------------------------------------------------------------------
--Inputs:
---rdfe: list of smoothness parameter r, number of interior two-faces numFaces, Euler characteristic e, and (optional) degree d
---edgeIncidences: a list recording the number of two-faces meeting at each interior edge
---intStarIncidences: for each interior vertex, the input numFaces, numslopeslist, required for vertexContribution
---bndryStarIncidences: for each boundary vertex, the input numFaces, numslopeslist, numpolytopes required for vertexContribution
-----------------------------------------------------------------------    
    r := rdfe_0;
    numFaces := rdfe_1;
    eulerC := rdfe_2;
    --number of interior vertices
    vint := length intStarIncidences;
    --main term coefficient
    tetC := eulerC+vint;
    --face term
    ntlist := apply(edgeIncidences,ed->min(r+2,ed));
    faceC := numFaces-sum(ntlist);
    --edge term appears below
    --vertex term
    ivertexT := sum apply(intStarIncidences,iv->vertexContribution(r,iv_0,iv_1,0));
    bvertexT := sum apply(bndryStarIncidences,bv->vertexContribution(r,bv_0,bv_1,bv_2));
    if length(rdfe)==3 then(
	S := QQ[d];
	--edge term
	edgeT := sum apply(ntlist, s-> edgeTerm(r,s,3,BaseRing=>S));
	return tetC*binomial(d+3,3)+faceC*binomial(d+2-r,3)+edgeT+ivertexT+bvertexT
	)else(
	d := rdfe_3;
	--edge term
	edgeT = sum apply(ntlist, s-> edgeTerm(r,d,s,3));
	return tetC*binom(d+3,3)+faceC*binom(d+2-r,3)+edgeT+ivertexT+bvertexT
	)
    )

trivariateLower(List,List) := (rd,F)->(
-----------------------------------------------------------------------
--Inputs:
---rd: list with r and d (or just r)
---F: list of tetrahedra (must be tetrahedral!)
-----------------------------------------------------------------------
    if all(F,f->#f==4) then(
	--get all two-faces
	allFaces := unique apply(flatten apply(F,f->subsets(f,3)),e->sort(e));
	--get interior dim 2 faces
	intFaces := select(allFaces, e->#select(F,f->isSubset(set e,set f))==2);
	--get number of interior dim 2 faces
	numFaces := length intFaces;
	--get boundary dim 2 faces
	bFaces := select(allFaces, e->#select(F,f->isSubset(set e,set f))==1);
	--get all edges
	allEdges := unique apply(flatten apply(F,f->subsets(f,2)),r->sort(r));
	--get boundary edges
	bEdges := unique apply(flatten apply(bFaces,s->subsets(s,2)),t->sort(t));
	--get interior edges
	intEdges := toList(set(allEdges)-set(bEdges));
	--number of interior edges
	numEdges := length intEdges;
	--get all vertices
	allVertices := sort unique flatten F;
	--get boundary vertices
	bVertices := sort unique flatten bFaces;
        --get interior vertices
	intVertices := toList(set(allVertices)-set(bVertices));
	--number of interior vertices
	numVertices := length intVertices;
	--get edgeIncidences list
	edgeIncidencesHash := hashTable apply(intEdges,e->{e,#select(intFaces,f->isSubset(e,f))});
	edgeIncidences := values edgeIncidencesHash;
	--get intStarIncidences
	intStarIncidences := apply(intVertices,v->(
		--get number of two-dim faces containing v
		nFv := #select(intFaces,f->member(v,f));
		--get edges containing v
		ev := select(intEdges,e->member(v,e));
		--get numslopes around each edge
		numslopesv := apply(ev,e->edgeIncidencesHash#e);
		{nFv,numslopesv}
		));
	--get bndryStarIncidences
	bndryStarIncidences := apply(bVertices,v->(
	    	--get number of two-dim faces containing v
		nFv := #select(intFaces,f->member(v,f));
		--get edges containing v
		ev := select(intEdges,e->member(v,e));
		--get numslopes around each edge
		numslopesv := apply(ev,e->edgeIncidencesHash#e);
		--get number of tetrahedra containing v
		numsimplicesv := #select(F,f->member(v,f));
		{nFv,numslopesv,numsimplicesv}
		));
	--Get Euler characteristic
	eulerC := length(F)-numFaces+numEdges-numVertices;
	--Plug data into previously defined trivariateLower
	if (length rd)==2 then(
	    rdfe := {rd_0,numFaces,eulerC,rd_1}
	    )else(
	    rdfe = {rd_0,numFaces,eulerC}
	    );
	return trivariateLower(rdfe,edgeIncidences,intStarIncidences,bndryStarIncidences)
	)else(
	return error "Must be simplicial for this form of input"
	)
    )
	
	
    