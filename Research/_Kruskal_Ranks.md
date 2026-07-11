# Maximal Kruskal ranks
## Kruskal ranks

Let $A \subset \mathbb{P}^n$ be a collection of points. The $d$th Kruskal rank $k_d(
A)$ of $A$ is the maximum value of $p$ so that every subset of $A$ of size $p$ imposes $p$ conditions on forms of degree $d$.

The $d$th Kruskal rank $k_d(A)$ is maximal if $k_d(A)=\min\{\binom{n+d}{n},|A|\}$.  To check that $k_d(A)$ is maximal, it suffices to check that every subset of $A$ of size $\binom{n+d}{n}$ imposes independent conditions on forms of degree $d$.  Equivalently, it suffices to check that $v_d(S)$ is linearly independent for every subset $S\subset A$, where $v_d$ is the $d$th Veronese map.

## The Veronese map

We first code the veronese map.
```
veronese = (U,d)->apply(compositions(#U,d),c->product apply(#U,i->U_i^(c_i)))
```

## Maximal Kruskal ranks

The running example in the final section of our paper is given by the columns of 
```
M=matrix {{2, 4, 0, 3, 5, 2, 1, 9, 4, 7, 6, 5, 5}, {9, 9, 1, 5, 3, 7, 2, 3, 5, 8, 2, 7, 1}, {4, 8, 6, 1, 0, 8, 6, 2, 8, 8, 4, 2, 2}, {0, 9, 1, 5, 5, 1, 6, 2, 3, 4, 0, 5, 7}}
```
To check that $k_1(A)$ is maximal, we check that every subset of $\binom{3+1}{3}=4$ columns of $M$ is linearly independent.
```
all(subsets(13,4),s->rank(M_s)==4)
```
To check that $k_2(A)$ is maximal, we check that $v_2(S)$ is linearly independent for every subset of size $\binom{3+2}{3}=10$ of $A$.
```
v2M=transpose matrix apply(entries transpose M,row->veronese(row,2))
```
```
all(subsets(13,10),s->rank(v2M)==10)
```
Finally, to check that $k_3(A)$ is maximal, we check that $v_3(A) is linearly independent.
```
v3M=transpose matrix apply(entries transpose M,row->veronese(row,3))
rank v3M
```
This verifies that $A$ is general in the sense of our paper.