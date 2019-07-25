--In this file we illustrate Chebotarev's Density Theorem
--Start with the polynomial you'd like to work with

--copy and paste this block of code into Macaulay2 to define 
--the function below.
--Input: a polynomial f over the integers and "numprimes", the number of primes
--to reduce f over.
--Output: a tally of how many times different factorization types of f occur
--among the first "numprimes"
tallyFactorizationStatistics:=(f,numprimes)->(
    cPrime := 1;
    counter := 0;
    factList := [];
    while counter<numprimes do(
	p := nextPrime(cPrime);
	F := (ZZ/p)[x];
	f = sub(f,F);
	T := factor f;
	ex := apply(#T,i->T#i#1);
	factList = append(factList,ex);
	cPrime = p+1;
	counter = counter+1;
	);
    tally factList
    )

--Here's how you can use the function.
--First, define the polynomial ring in one variable over the integers:
R=ZZ[x];
--next, define your favorite polynomial:
f1=x^4+5*x^3+x^2-2;
f2=x^4+3*x^2-1;
f3=x^4+3*x^2+1;
f4=x^5+10*x^3-10*x^2+35*x-18;
f5=x^5-3*x-1;
--now, use the function defined above.
tallyFactorizationStatistics(f1,10000)
