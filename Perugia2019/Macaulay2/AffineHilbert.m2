--A script to compute the affine Hilbert function
affineHilbertFunction=(d,I)->(
    S := ring(I);
    afH := 0;
    for i from 0 to d do(
	afH = afH+numcols basis(i,S/I)
	);
    afH
    )

--A script to display affine Hilbert function up to a certain point
affineHilbertFunctionTable=(d,I)->(
    S := ring(I);
    L := apply(d+1,i->numcols basis(i,S/I));
    A := accumulate(plus,prepend(0,L));
    tp := prepend("d",toList(0..d));
    bm := prepend("AffineHilbert(d)",A);
    netList {tp,bm}
    )