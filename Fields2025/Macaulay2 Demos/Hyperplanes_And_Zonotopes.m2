# Hyperplanes and Zonotopes
## Loading packages

For this tutorial we'll need to load both the "HyperplaneArrangements" package and the "AlgebraicSplines" package.

```
loadPackage("HyperplaneArrangements",Reload=>true);
loadPackage("AlgebraicSplines",Reload=>true);
```

The "HyperplaneArrangements" package loads the package "Polyhedra", which we will also need.

## Hyperplane arrangements

A hyperplane arrangement $A\subset\mathbb{R}^n$ is a union
\[
A=\bigcup_{i=1}^k H_i
\]
where $H_1,\ldots,H_k$ are hyperplanes -- i.e. $H_1,\ldots,H_k$ are each the zero locus of a linear form $\alpha_1,\ldots,\alpha_k$.  An easy way to refer to a hyperplane arrangement is to refer to the list $\{\alpha_1,\ldots,\alpha_k\}$ (this is the input we will use for the hyperplane arrangements package) or as the product
\[
Q(A)=\prod_{i=1}^k \alpha_i,
\]
which is the *defining equation* of $A$.

The examples of splines on quadrants and octants come from *boolean* arrangements -- these are the arrangements with defining equations $xy$ in $\mathbb{R}[x,y]$ and $xyz$ in $\mathbb{R}[x,y,z]$.  We can define these as follows (note that you first need to define the ambient ring).

```
S1=QQ[x,y];
A1=arrangement {x,y}
```

```
S2=QQ[x,y,z];
A2=arrangement {x,y,z}
```

One of the things that is useful to know is where the hyperplanes intersect each other.  In three dimensions, since we will be concerned with *central* hyperplane arrangments (all hyperplanes pass through the origin), we really only need to know about intersections along lines.  We get this information with the command *flats*

```
--flats(d,A) returns the all sets of hyperplanes that intersect in a codimension d linear space
flats(0,A1)--by convention, the intersection of no hyperplanes corresponds to the ambiend space
flats(1,A1)--the codimension one `intersections' are just the hyperplanes
flats(2,A1)--there is one codimension two intersection - the origin
```
```
flats(0,A2)
flats(1,A2)
flats(2,A2)
flats(3,A2)
```
These are more readable if converted into lists.  For instance:
```
F2=flats(2,A2);
apply(F2,f->toList f)
```


## Chambers of hyperplane arrangements
To define splines, we need some way to access and record the *chambers* of a hyperplane arrangement $A\subset\mathbb{R}^n$ (we'll only be looking at $n=2,3$).  A chamber of a hyperplane arrangement is the closure (in the Euclidean topology) of a connected component of $\mathbb{R}^n\setminus A$.  For instance, the chambers of the boolean arrangement $xy$ in $\mathbb{R}^2$ are the four quadrants, and the chambers of $xyz$ in $\mathbb{R}^3$ are the eight octants.

The chambers of a hyperplane arrangement are the $n$-dimensional cones of something called a *fan*.  We'll denote the fan associated to an arrangement $A$ by $\Sigma^A$.

What we need to define splines on the chambers (equivalently, on the fan $\Sigma^A$), is the data of the chambers along with which chambers intersect along a codimension one linear space (i.e. *adjacent* chambers).  That is, what we need is a *dual graph*.

Luckily, there is such an object!  In fact, there is a whole polytope, called a *zonotope*, associated to the hyperplane arrangement $A$, which we will denote by $Z(A)$.  The vertices of $Z(A)$ correspond to the $n$-dimensional cones of $\Sigma^A$ (the chambers), the edges of $Z(A)$ correspond to pairs of adjacent chambers (i.e. cones of codimension one in $\Sigma^A$), and in general the $i$-dimensional faces of $Z(A)$ correspond to $(n-i)$-dimensional cones of $\Sigma^A$.

Also fortunately, $Z(A)$ is not so hard to define.  Suppose that $A\subset\mathbb{R}^n$ is defined by linear forms $\alpha_1,\ldots,\alpha_k$.  Let $v_1,\ldots,v_k\in\mathbb{R}^n$ be the normal vectors of $\alpha_1,\ldots,\alpha_k$ (that is, these record the coefficients of the linear forms).  Let $[-v_i,v_i]$ be the line segment from $-v_i$ to $v_i$ for $i=1,\ldots,k$.  Then $Z(A)$ is the *Minkowski sum*
\[
Z(A)=[-v_1,v_1]+[-v_2,v_2]+\cdots+[-v_k,v_k]
\]
of the $k$ line segments.

Example.  Suppose $Q(A)=xy$.  Then $v_1=(1,0)$ and $v_2=(0,1)$.
```
v1=transpose matrix{{1,0}}
M1=(-v1)|v1
L1=convexHull(M1)
```
The last command (from the 'Polyhedra' package) takes the convex hull of the columns of M1.
```
vertices L1 --vertices of the polytope (line segment!) L1
```
Since L1 is a polytope, it is defined by inequalities, which we can get as follows.
```
facets L1 --this records the inequalities that define L1 as a polytope
```
Now do the same for v2
```
v2=transpose matrix{{0,1}}
M2=(-v2)|v2
L2=convexHull(M2)
```
Minkowski sum is just ordinary sum!
```
Z=L1+L2
```
Let's call up the vertices of Z.
```
vertices Z --four vertices, one for each chamber of $xy$
```
And the facets.
```
facets Z --four inequalities, one for each ray of the fan \Sigma^A
```

All these relationships between hyperplane arrangements and zonotopes (and a great deal more) can be found in Ziegler's excellent book *Lectures on Polytopes*, chapter seven.

## Generating the input for the generalized splines function

In order to compute splines on the fan $\Sigma^A$ induced by a hyperplane arrangement $A\subset \mathbb{R}^n$ with a smoothness distribution $r:\Sigma^A_{n-1}\rightarrow \mathbb{Z}$, we need to extract the data of the vertices and edges of the zonotope $Z(A)$, along with the linear forms corresponding to each edge of the zonotope.  We can get this data as follows: if $v,w$ are adjacent vertices of $Z(A)$, then $v-w=\pm 2v_i$ for some vector $v_i$ orthogonal to a hyperplane of $A$.  Click on the following code to define the function "zonotopeEdgeLabels" that puts this all together.  We'll see how to use this function on the next screen.  All functions defined in this tutorial can also be found in the file "Scripts_for_splines_on_hyperplanes.m2", which can be loaded independently to the web interface.

```
--Inputs: A: a hyperplane arrangement A
--Outputs: A sequence (E,I) where
--        E records the edges between vertices of Z(A)
--        I records the edge labels using the smoothness distribution r

zonotopeEdgeLabels = method()
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
    while vtxcnt<(#V)-1 do(
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
```

## Using the zonotopeEdgeLabels function

Let's recover the splines on the boolean arrangement $xy$ and $xyz$ using the zonotopeEdgeLabels function.

```
S=QQ[x,y]
A1=arrangement{x,y}
(E1,L1,Z1)=zonotopeEdgeLabels(A1)
Spl1=generalizedSplines(E1,L1)--this returns the continous splines (no higher derivatives)
```
Z1 in the code above is the actual zonotope.

```
S=QQ[x,y,z]
A2=arrangement{x,y,z}
(E2,L2,Z2)=zonotopeEdgeLabels(A2)--this returns the continous splines (no higher derivatives)
Spl2=generalizedSplines(E2,L2)
```

## Changing the smoothness distribution
In our project, we will be changing the smoothness distribution across hyperplanes.  It will be useful to have a short script which goes through the list of edge labels and changes the exponents according to the smoothness distribution on the hyperplanes.  The following function does this.

```
assignSmoothness = method()
assignSmoothness = (A,r,I)->(
    hyps := A#(first keys A);--get the list of hyperplanes
    smthHash := hashTable apply(length r,i->{hyps_i,r_i});--create a hash table that associates each linear form of the hyperplane to the desired smoothness
    apply(I,h->h^(smthHash#h+1))--cycle through I, raising each linear form to the appropriate power
    )
```

Here it is in operation.

```
S=QQ[x,y,z]
A=arrangement {x,y,z}
(E,L,Z)=zonotopeEdgeLabels(A)
r={3,4,5}
I=assignSmoothness(A,r,L)
```
Now we can feed $E$ and $I$ into generalizedSplines.

```
generalizedSplines(E,I)
```
