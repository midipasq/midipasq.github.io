--Implementation of the Borel sort algorithm in "Koszul multi-Rees algebras of principal L-Borel ideals" 
--Any references to lemmas or theorems refer to this paper

--Input: a monomial m=x_1^a_1...x_n^a_n and an integer i, 0<i<n+1
--Output: cumulative exponent a_(n-i)+a_(n-i+1)+...+a_n
cexp=(i,m)->(
    S := ring(m);
    n := length gens S;
    expts := flatten exponents(m);
    inds := apply(n-i,j->i+j);
    sum apply(inds,i->expts_i)
    )

--Input: Two monomials M,N.
--Output: True if M can be obtained from N by borel moves, otherwise false
borelOrder=(M,N)->(
    expM := flatten exponents(M);
    expN := flatten exponents(N);
    dm := sum expM;
    dn := sum expN;
    bool := false;
    --If M and N have different degrees then output 'false'
    if dm==dn then(
	bool = true;
	cntr := 1;
	n := length expN;
	--check cumulative exponents
	while (bool and cntr<n) do(
	    ind := apply(cntr,i->n-i-1);
	    sM := sum apply(ind, j-> expM_j);
	    sN := sum apply(ind, j-> expN_j);
	    if sM>sN then bool=false;
	    cntr = cntr+1
	    ));
    bool
    )

--Input: A monomial M, variable v (not the largest), integer q.
--Output: least monomial under Borel order in Borel(M) not divisible by variables 'less than' v and having exponent q on v (see Lemma 4.2)

borelMins=(M,v,q)->(
    S := ring(M);
    V := gens S;
    s := position(V,i->i==v);
    --return error if there is no monomial in Borel(M) with power of q on v
    if cexp(s,M)<q then(
	return("No monomial in Borel("|toString(M)|") is divisible by "|toString(v^q))
	);
    if v==(first gens S) then(
	return("Input variable should not be the largest variable.")
	)else(
	Mexp := flatten exponents M;
	Nexp := Mexp_(apply(s-1,i->i));
	N := (product apply(s-1,i->V_i^(Nexp_i)))*v^q;
	vs := V_(s-1);
	D := (first flatten degree M)-(first flatten degree N);
	N = N*vs^D;
	);
    N
    )

--Input: A target multi-degree MD (\mu in the paper) and a monomial M
--Output: the borel sort of MD into k factors from Borel(M) (Algorithm 1)
borelSort=(M,MD)->(
    --get degree of M
    d := first degree M;
    --get k so that MD is in Borel(M^k)
    k := floor((first degree MD)/d);
    --this error will only ever be called on the first level of recursion
    if not borelOrder(MD,M^k) then(
	return "No factorizations possible"
	);
    --get internal copies of the ambient ring of MD and M
    S := ring(MD);
    MD = sub(MD,S);
    M = sub(M,S);
    --get the variables of S
    VS := gens(S);
    --the following if clause is triggered when M has degree 1 (should be un-necessary)
    --if d==1 then(
    --return(flatten apply(length VS,i->apply((flatten exponents MD)_i,j->VS_i)))
    --);
    --get the variable of largest index dividing MD
    lvar := last(select(VS,v->(MD%v==0)));
    --get the exponents of MD
    expMD := flatten exponents MD;
    --get the power of the variable lvar in MD
    A := last(select(expMD,i->i!=0));
    --The following `if' clause is the stopping criterion corresponding to the first if clause in Algorithm 1
    if MD==lvar^(k*d) then(
	--Update mu_1,...,mu_k, stored in the list bSort
	bSort := apply(k,i->lvar^d)
	)else(
	q := floor(A/k);
	r := A%k;
	if r>0 then(
	    --Mup is the monomial M_up in Algorithm 1
	    Mup := sub((1/lvar^q)*borelMins(M,lvar,q),S);
	    --MDup is the monomial \mu_up in Algorithm 1
	    MDup := last select(flatten entries borel matrix{{Mup^(k-r)}}, b->MD%b==0);
	    --Mdown is the monomial M_down in Algorithm 1
	    Mdown := sub((1/lvar^(q+1))*borelMins(M,lvar,q+1),S);
	    --MDdown is the monomial \mu_down in Algorithm 1
	    MDdown := sub(MD/(MDup*lvar^A),S);
	    --Next line is Borel Sort of M_up, \mu_up
	    bsUpper := borelSort(Mup,MDup);
	    --Next line is Borel sort of M_down, \mu_down
	    if d == (q+1) then(
		bsLower := apply(r,i->1)
		)else(
		bsLower = borelSort(Mdown,MDdown)
		);
	    --Next line updates mu_1,...,mu_k, stored in the list bSort
	    bSort = apply(k,i->(if i<k-r then lvar^q*bsUpper_i else lvar^(q+1)*bsLower_(i-k+r)))
	    )else(
	    --Mleft is the monomial M_left in Algorithm 1
	    Mleft := sub((1/lvar^q)*borelMins(M,lvar,q),S);
	    --MDleft is the monomial \mu_left in Algorithm 1
	    MDleft := sub(MD/lvar^A,S);
	    --Next line is Borel sort of Mleft, \mu_left
	    bsLeft := borelSort(Mleft,MDleft);
	    --Next line updates mu_1,...,mu_k, stored in the list bSort
	    bSort = apply(k,i->lvar^q*bsLeft_i)
	    )
	);
    --return the BorelSort
    bSort
    )
