--The following code defines the function 'cat', which produces the catalecticant matrix of a degree.  Copy and paste into a macaulay2 session to define this function.
--Input: a form f and a degree k.
--Output: a matrix for the map from S_k to S_(d-k) given by letting polynomials in S_k act by partial differentiation on f
cat=(f,k)->(
    S := ring(f);
    DSk := basis(k,S);
    d := first degree f;
    MonsSource := flatten entries DSk;
    MonsTarget := flatten entries basis(d-k,S);
    matrix apply(MonsSource,B->flatten entries transpose ((coefficients(diff(B,f),Monomials=>MonsTarget))_1))
    )

--Define ring containing necessary variables.
RT=QQ[A..V,a..f]
RF=frac(QQ[A..V,a..f])[x,y,z];
--The defining equation of an irreducible multi-arrangement in P^2, in standard form
HA=x*y*z*(x+y+z)*(a*x+b*y+c*z)*(d*x+e*y+f*z)
--The catalecticant matrix of this form (first need to copy and paste function above)
CAT=sub(cat(HA,3),RT)

--Force rank 7 or lower by multiplying by rank 3 matrix
use RT
M=((id_(RT^3))||matrix({{A,B,C},{D,E,F},{G,H,V},{J,K,L},{M,N,O},{P,Q,R},{S,T,U}}))
INLIST=matrix({flatten entries(CAT*M)})

--Square this system by taking 25 random linear combinations of the equations in INLIST
RAN=random(RT^30,RT^25)
EQLIST=flatten entries(INLIST*RAN)

--format input for Bertini
EQLIST1=apply(length EQLIST,i->concatenate("F",toString(i+1),"=",toString(EQLIST_i),";"));
SLIST="";
for i from 0 to (length EQLIST1-1) do(
    if i==0 then(
	SLIST=concatenate(SLIST,EQLIST1_i)
	)else(
	SLIST=concatenate(SLIST,concatenate("\n",EQLIST1_i))
	))
--After running the above code to define SLIST, copy the string SLIST and paste it into a text file to run Bertini on.
--A sample of such a file is in the file input_cubic_annihilator.
SLIST
