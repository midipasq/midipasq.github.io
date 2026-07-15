# Verifying that a given list of polynomials defines an aCM sextic of genus three

We verify that the following ideal defines a smooth aCM sextic curve in \(\mathbb{P}^4\) of degree $6$ and genus $3$.

```
R=QQ[w,x,y,z];
M=matrix({{w, x, y, z},{x+2*y+3*z, 2*w+7*y+9*z, 21*x+22*y, 23*y+17*z},{13*x+2*y-3*z, w-4*x+y-13*z, w+2*x-y+5*z, 3*w+2*x-7*y}})
I=minors(3,M)
```

First we verify that $I$ defines a curve.

```
isPrime(I)
```

```
codim(I) --should be two to define a curve in P^3
```

We now verify that the degree and genus are correct.

```
genus(I) --should be three
```

```
degree(I) --should be six
```

We verify that I defines a smooth variety.

```
isSmooth(I)
```

Finally, we check that $I$ defines an aCM curve.

```
pdim(coker gens I)==codim(I) --should be true
```