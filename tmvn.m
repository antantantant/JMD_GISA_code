function z=tmvn(bd)

phi=normcdf(bd(2));
plo=normcdf(bd(1));
r=rand();
r=plo+(phi-plo)*r;
z=norminv(r);