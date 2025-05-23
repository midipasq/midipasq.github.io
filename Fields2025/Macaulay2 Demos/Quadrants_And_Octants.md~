{{ page.markdown }}

# First examples: Splines on quadrants and octants
## Setup: Loading the package and generalizedSplines function

First we will need to load the AlgebraicSplines package.  Click the 'Start' button above the right-hand panel to start a Macaulay2 session, and then click the link below.

```
loadPackage("AlgebraicSplines",Reload=>true)
```

The basic function we will use is generalizedSplines.  We can read the documentation on this function as follows:

```
help generalizedSplines
```

The basic input for generalizedSplines is a graph recorded as a list of edges (the vertices are assumed to be integers).

## Splines on the quadrants in $\mathbb{R}^2$
Let's start by computing the module of splines on the fan in $\mathbb{R}^2$ whose four maximal cones are the four quadrants.

We label the quadrants as: first, the upper right, second, the upper left, third, the lower right, and fourth, the lower left.  With this labeling, the adjacent pairs of quadrants are shown by the edge list E defined below. (Notice that lists are defined using curly braces in Macaulay).

```
E={{1,2},{1,3},{2,4},{3,4}}
```

If we don't want the output for E to get printed, we can put a semicolon at the end:

```
E={{1,2},{1,3},{2,4},{3,4}};
```

Now let's give the edge labels.  We need to start by defining an ambient polynomial ring.  We'll use the polynomial ring in two variables over the rationals.

```
R=QQ[x,y]
```

Now we'll do the edge labels (start with the smoothness distribution of 0 on all edges):

```
I={x,y,y,x}
```

How did we get this?  The quadrants 1 and 2 meet along the y-axis, where x vanishes, so x is the label we give to the edge {1,2}, and so on.

Now we've defined the input we need for the generalizedSplines function, so we can go ahead and compute.

```
M=generalizedSplines(E,I)
```

Macaulay presents the result as an image of a matrix.  We can get the generators as a matrix like this:
```
gens M
```

Taking the determinant gives

```
det gens M
```

Notice that, up to constant multiple, this is the product of the edge labels.  We can change the smoothness distribution as follows.

```
I={x^3,y^2,y^2,x^3}--corresponds to a smoothness distribution of 2,1,1,2 across the rays
M=generalizedSplines(E,I)
```
Macaulay can recognize this is a free module as follows
```
prune M
```
See if you can match up the degree information with the columns of the matrix.

## Splines on the octants in $\mathbb{R}^3$
Now let's do a three-dimensional example where the fan consists of the eight octants in $\mathbb{R}^3$.

Let's use the convention that +++ refers to the octant where all coordinates are non-negative and +-+ refers to the octant where the $x$ and $z$ coordinates are non-negative while the $y$ coordinate is non-positive, and so on.

We list the octants in order as +++,-++,+-+,++-,--+,-+-,+--,---.  Now we list pairs of adjacent octants using this ordering.
```
E={{1,2},{1,3},{1,4},{2,5},{2,6},{3,5},{3,7},{4,6},{4,7},{5,8},{6,8},{7,8}}
```
The polynomial ring is now in three variables.
```
R=QQ[x,y,z]
```

Let's again do just the continuous smoothness distribution (assign $0$ to all edges).  We can get the right linear form by looking at which sign changes between the adjacent octants.
```
I={x,y,z,y,z,x,z,x,y,z,y,x}
```
Let's compute the generalized splines.
```
M=generalizedSplines(E,I);
```
Get the matrix of generators:
```
gens M
```
Observe the determinant is again the product of the edge labels.
```
det gens M
```
So this is a free spline module, which we can also recognize by using the prune command again.
```
prune M
```
