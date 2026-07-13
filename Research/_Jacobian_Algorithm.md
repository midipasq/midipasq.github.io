# Jacobian Algorithm
## Idea of the algorithm

Let $A \subset \mathbb{P}^3$ be a collection of points in general position (see the tutorial on Kruskal ranks) and $S=\mathbb{K}[x_0,\ldots,x_3]$ the homogeneous coordinate ring of $\mathbb{P}^n$.  The Jacobian algorithm, which is Algorithm 3 in Section 6 of our paper, gives a criterion to verify that $A$ does not lie on a curve $C$ of a specific degree $d$ and genus $g$.

Let $I_C$ be the ideal of the curve.  The main observation to get the algorithm started is that, since we know the degree and genus $C$, then we know the Hilbert polynomial of the curve is $HP_{S/I_C}(n)=dn+1-g$.  Moreover,by the Gruson-Lazarsfeld-Peskine theorem, \(\mathrm{reg}(S/I_C)\le \deg(C)-\mathrm{codim}(C)\).  Hence
\[
\mathrm{dim}(S/I_C)_n=dn+1-g
\]
for \(n\ge d-2\).  Moreover, if \(S/I_{C}\) is Cohen-Macaulay, \(\mathrm{dim}(S/I_C)_n=dn+1-g\) for \(n\ge d-3\).

Let $n\ge d-2$ (or $n\ge d-3$ if $S/I_C$ is Cohen-Macaulay) and assume that $|A|\le \binom{n+3}{3}$.  Since we have the coordinates of the points in $A$, we also have access to the ideal $I_A$ of the set of points.  Let
\[
e=\mathrm{dim} (I_C)_n=\binom{n+3}{3}-dn+g-1
\]
and
\[
s=\mathrm{dim}(I_A)_n=\binom{n+3}{3}-|A|
\]
(which follows since $A$ is in general position).  The idea of the Jacobian algorithm is to scan all $e$-dimensional subspaces of $(I_A)_n$ to see if any locally cut out a 1-dimensional variety at each of the points of $A$.  If not, then there is no integral curve C of degree $d$ and genus $g$ that contains $A$.

## Illustrating the algorithm for thirteen points: the ideal of the points

We continue with the running example in Section 6 of our paper, which is the set of points $A$ whose coordinates are given by the columns of the matrix $M$ below:

```
K=QQ; --choose a field to work over
S=K[x_0,x_1,x_2,x_3];
M=matrix {{2, 4, 0, 3, 5, 2, 1, 9, 4, 7, 6, 5, 5}, {9, 9, 1, 5, 3, 7, 2, 3, 5, 8, 2, 7, 1}, {4, 8, 6, 1, 0, 8, 6, 2, 8, 8, 4, 2, 2}, {0, 9, 1, 5, 5, 1, 6, 2, 3, 4, 0, 5, 7}};
```

We again need the ideal of $A$, so we create the script to compute that ideal

```
--Inputs:
--S: the ring to contain the ideal of points, which should have the same number of variables as the coordinates of the points
--M: a matrix whose rows are coordinates of the points

pointIdeal = (S,M)->intersect apply(entries M,r->minors(2,matrix {gens S,r}));
I = pointIdeal(S,transpose M);
```

We look for an aCM integral sextic curve $C$ of genus three which contains $A$.  By the Gruson-Lazarsfeld-Peskine bound, the regularity of $S/I_C$ is at most $3$.  Thus, since $S/I_C$ is Cohen-Macaulay,
\[
\dim(S/I_C)_n=6d-2
\]
for $d\ge 3$.  So we take $n=3$ (in the description on the first slide) and we need a basis for the cubics which vanish on $A$.

```
CUB = super basis(3,I); --CUB is a basis for the cubics that vanish on A
```

The cubics that vanish on $C$ form a vector space of dimension $\binom{3+3}{3}-(6\cdot 3-2)=4$.  Our remaining task is to scan through the four-dimensional subspaces of $(I_C)_3$ to check if any such subspace vanishes on a curve.

## Illustrating the algorithm for thirteen points: the first patch

The four-dimensional subspaces of the seven dimensional space of cubics in \(I_A\) form a Grassmannian of type \(\mathbb{G}(4,7)\).  To see if there is any four-dimensional subspace of the space of cubics of \(I_A\) that could vanish on a curve (which is necessary if $A$ is contained in an aCM integral sextic curve of genus three), we run through the affine charts of $\mathbb{G}(4,7)$.  We start by illustrating the procedure for a single affine chart.  The following code produces a $7\times 4$ matrix $B$ which we will multiply by the basis of $(I_A)_3$.  `Section 6' refers to section 6 of our paper.

```
PAR=K[(symbol b)_(1,1)..(symbol b)_(3,4)];--polynomial ring with 12 variables to fill a three by four submatrix
ID=id_(PAR^4);
Ap=(matrix table(toList(1..3),toList(1..4),(i,j)->b_(i,j)));
B=ID||Ap--the 'B' matrix in section 6 which we use to parameterize a four-dimensional subspace
```

Next we define a polynomial ring that contains all the variables and form a list of substitutions for the x-variables corresponding to the coordinates of the points.

```
BRI=PAR[(symbol x)_0..(symbol x)_3];
ePTS=entries sub(transpose M,BRI);--substitute points into the ring BRI
SUBS=apply(ePTS,p->apply(length p,i->BRI_i=>p_i));--get list of 13 coordinate substitutions to apply to the Jacobian matrix
```

Next take the Jacobian of the matrix of cubics CUB which span the degree three part of the ideal $I_A$.  Apply the thirteen point substitutions and multiply by the matrix B to get the Jacobian of a general element of a four-dimensional span of the cubic generators of $I_A$, on the given coordinate chart.

```
F=sub(CUB,BRI);--the 'F' matrix in section 6.
JACF=jacobian F; --The Jacobian of the cubic generators of $I_A$.
JACLISTF=apply(SUBS,s->substitute(substitute(JACF,s),PAR)); --Evaluate the Jacobian at the thirteen points
JACLIST=apply(JACLISTF,mn->mn*B); --multiply by B to get the Jacobian of the general element of a four-dimensional span of the cubic generators of $I_A$, evaluated at each point;
```

If there is an aCM sextic curve $C$ which vanishes at each point of $A$, then there will be a choice of $b_{1,1},\ldots,b_{3,4}$ so that the rank of each matrix in JACLIST is two; that is, the $3\times 3$ minors of each matrix in JACLIST vanish.  We next take the ideal generated by all $3\times 3$ minors of every matrix in JACLIST.

```
J=sum apply(13,i->minors(3,JACLIST_i));
```

According to the Nullstellensatz, if $1$ is in this ideal, there is no choice of $b_{1,1},\ldots,b_{3,4}$ so that the rank of each matrix in JACLIST is two.  That is, there is no aCM sextic curve $C$ of genus three which contains $A$.

```
1%J
```

Since $1$ is in the ideal $J$, there is no aCM sextic curve $C$ of genus three which contains $A$ (on the affine patch that we set up).

## Illustrating the algorithm for thirteen points: repeat for all affine patches

Finally, we run the same code over all affine patches to ensure that we did not miss anything.  This may take half a minute or so.

```
patches = apply(subsets(7,4),s->(
	cs = toList(0..6)-set(s);
	basischoice = F_s|F_cs;
	JACF=jacobian basischoice; --The Jacobian of the cubic generators of $I_A$.
	JACLISTF=apply(SUBS,s->substitute(substitute(JACF,s),PAR)); --Evaluate the Jacobian at the thirteen points
	JACLIST=apply(JACLISTF,mn->mn*B);
	--throw all jacobian conditions into an ideal
	J=sum apply(13,i->minors(3,JACLIST_i));
	1%J
	));

all(patches,i->i==0_PAR)
```

This shows that there is no four-dimensional subspace of $(I_A)_3$ which locally defines a curve at each of the points of $A$.  Thus there is no aCM sextic of genus $3$ containing $A$.