--code to produce derivation complex
loadPackage "HyperplaneArrangements"

--getting relations--
matRels=method()
------
--Input: M, a list of rows of a matrix, and R, a choice of rows of M
--Output: N, a list of relations on the R rows of M,
------encoded in the space of relations on all rows
------
matRels(List,List):=List=>(M,R)->(
    rs := length M;
    cs := length (M_0);
    c := rs-(#R);
    if c==0 then(
	return entries transpose gens ker transpose matrix M
	)else(
	MN := transpose matrix apply(rs,i->(
		if member(i,R) then (
		    M_i
		    )else(
		    apply(cs,j->0)
		    )));
	S :=ring(M_0_0);
	Id :=id_(S^c);
	counter := 0;
	i := 0;
	Ad := {};
	while i<rs do(
	    if member(i,R) then(
		Ad = append(Ad,apply(c,j->0));
		i = i+1
		) else (
		Ad = append(Ad,flatten entries (Id_counter));
		counter = counter+1;
		i = i+1
		));
	return entries transpose gens ker(MN||(transpose matrix Ad))
	)
    )


--the main function for formality--
formalityComplex=method()

------
--Input: A, an arrangement
--Output: A sequence (C,H)
-----C: Chain complex so that if A is k-formal then H^i(C)=0 for i=1,..,k-1
-----H: List of hash tables encoding flats which correspond to rows of the ith differential,
---------along with what rows of the matrix correspond to that flat
------

formalityComplex(Arrangement):=Sequence=>A->(
    S := A.ring;
    r := rank A;
    Fl := flats A;
    --initialize the list that will become a hash table 
    --to lookup what rows of the matrices belong to which flats
    M1 := entries transpose coefficients A;
    H1 := new HashTable from prepend(codim => 1, apply(length(Fl_1),i->((Fl_1)_i=>{i})));
    relInfo := {H1};
    --initialize the list that will become the maps in the
    --chain complex for formality
    Mats := {map(S^(length M1),S^(length (M1_0)),matrix M1)};
    --loop creating relInfo and Mats
    cod := 2;
    currentRank := 1;
    while currentRank>0 and cod<=r do(
	F := Fl_cod;
	MP := Mats_(cod-2);
	HP := relInfo_(cod-2);
	fp := select(keys HP,i->((class i))===Flat);
	ind := 0;
	flrows := 0;
	newMat := {};
	newH := {};
	while ind<length F do(
	    cFlat := F_ind;
	    covering := select(fp,g->isSubset(g.flat,cFlat.flat));
	    rows := join(apply(covering, c->HP#c));
	    newM := matRels(entries MP,flatten rows);
	    p := length(newM);
	    --fork depending on the existence of relations--
	    if p==0 then(
		ind = ind+1;
		)else(
		--indices of rows corresponding to cFlat
		newH = append(newH,cFlat=>toList(flrows..(flrows+p-1)));
		--join relations to newMat
		newMat = join(newMat,newM);
		flrows = flrows+p;
		ind = ind+1;
		)
	    );
	--fork depending on existence of relations
	if length(newMat)==0 then(
	    currentRank = 0;
	    )else(
	    --update relinfo and Mats
	    relInfo = append(relInfo,new HashTable from prepend(codim=>cod,newH));
	    Mats = append(Mats,map(S^(length newMat),S^(length (newMat_0)),matrix newMat));
	    cod = cod+1;
	    )
	);
    --create chain complex--
    C := chainComplex(reverse Mats);
    --return shifted chain complex so taking cohomology works well, along with relInfo
    (C[cod-1],relInfo)
    )

--creating the kernel which will be fed to derivationComplex--
kernelDerivationComplex=method()

-----------
--Input: A, an arrangement, m, a list of multiplicities on hyperplanes, and Q, output of formalityComplex(A)
--Output: J, a sub-chain complex of formalityComplex(A)
-----------

kernelDerivationComplex(Arrangement, List, Sequence):=ChainComplex=>(A,m,Q)->(
    S := A.ring;
    ell := numgens S;
    (C,H) := Q;
    Hyps := A#((keys A)_0);
    --intialization--
    kerns := {image transpose matrix{apply(ell,i->0_S)},directSum apply(length Hyps,i->image matrix{{(Hyps_i)^(m_i)}})};
    --while loop to inductively construct kernel modules--
    current := 2;
    while current<= (length C) do(
	Jprev := kerns_(current-1);
	--current relation info
	cRI := first select(H,h->(h#codim==current));
	--current flats
	ciFlats := select(keys cRI,i->((class i))===Flat);
	--sort them properly
	ordRows := sort(apply(ciFlats,F->(cRI#F)));
	cFlats := apply(ordRows,i->first select(ciFlats,F->(cRI#F==i)));
	--previous differential--
	pdiff := (C.dd)_(-current+1);
	newKerns := directSum apply(cFlats,F->(
		image((pdiff^(cRI#F)*(gens Jprev)))
	    ));
	kerns = append(kerns,newKerns);
	current = current +1;
	);
    Maps := apply(length C,i->inducedMap(kerns_(i+1),kerns_i,C.dd_(-i)));
    CJ := chainComplex(reverse Maps);
    CJ[length C]
    )
    

-----------
--Input: A, an arrangement, and m, a list of multiplicities on hyperplanes
--Output: J, a sub-chain complex of formalityComplex(A)
-----------

kernelDerivationComplex(Arrangement, List):=ChainComplex=>(A,m)->(
    Q := formalityComplex(A);
    kernelDerivationComplex(A,m,Q)
    )

kernelDerivationComplex(Arrangement):=ChainComplex=>A->(
    n := #flats(1,A);
    kernelDerivationComplex(A,apply(n,i->1))
    )
    
--the derivation complex
derivationComplex=method()

-------
--Input: A, an arrangement, and m, a multplicity
--Output: A chain complex whose homologies encode freeness of (A,m)
-------

derivationComplex(Arrangement,List):=ChainComplex=>(A,m)->(
    (C,H) := formalityComplex(A);
    CJ := kernelDerivationComplex(A,m,(C,H));
    coker(inducedMap(C,CJ))
    )

derivationComplex(Arrangement):=ChainComplex=>A->(
    n := #flats(1,A);
    derivationComplex(A,apply(n,i->1))
    )