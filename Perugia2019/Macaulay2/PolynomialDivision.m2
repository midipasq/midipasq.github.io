--Division algorithm from CLO p.65 implemented in Macaulay2
--Inputs: polynomial f and list L of polynomials to divide f by
--Outputs: quotient list quot and remainder rem of division.  Satisfies sum q_iL_i+rem=f,
--and none of the terms of rem are divisible by the leading terms of polynomials in L
--Copy and paste the following code into Macaulay2:
divisionAlgo=(f,L)->(
    --get the ring in which polynomials are defined
    S := ring f;
    --initialize quotient list and remainder to zero
    quot := new MutableList from apply(L,i->0);
    rem := 0;
    p := f;
    while p != 0 do(
	i := 0;
	divisionoccurred := false;
	while i<=(length(L)-1) and not divisionoccurred do(
	    diffi := (flatten exponents leadTerm(p))-(flatten exponents leadTerm(L_i));
	    if all(diffi,i->i>=0) then(
		nqt := sub(leadTerm(p)/leadTerm(L_i),S);
		quot#i = quot#i+nqt;
		p = p-nqt*L_i;
		divisionoccurred = true;
		)else(
		i = i+1;
		));
	if not divisionoccurred then(
	    rem = rem+leadTerm(p);
	    p = p-leadTerm(p);
	    );
	);
    (toList quot,rem)	
    )

----------------------------------
--Illustration:-------------------
----------------------------------
--copy and paste the above function into
--Macaulay2, then execute the following 
--lines of code.  Change the term order
--and see what kind of results you get.
S=QQ[x,y,MonomialOrder=>GLex]--To switch to lex order use MonomialOrder=>Lex
f=x^7*y^2+x^3*y^2-y+1
L={x*y^2-x,x-y^3}
(q,r)=divisionAlgo(f,L)--this assigns the quotient list to the variable q and the remainder to the variable r
(q,r)=divisionAlgo(f,reverse L)
--Check:this should evaluate to f
(sum apply(length q,i->q_i*L_i))+r

f=x^9*y^4-3*x^3*y^2+5*x*y^4
L={x*y^2-x,x-y^3}
(q,r)=divisionAlgo(f,reverse L)
----------------------------------
----------------------------------
----------------------------------


--The function quotientRemainder actually computes a Grobner basis in the background
--and uses that to minimize the remainder
S=QQ[x,y,MonomialOrder=>GLex]
f=matrix{{x^7*y^2+x^3*y^2-y+1}}
M=matrix{{x*y^2-x,x-y^3}}
(Q,R)=quotientRemainder(f,M)
--Check:
f==M*Q+R

