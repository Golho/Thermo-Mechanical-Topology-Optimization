%    This is the file beammain.m
%    which is used as a main program for
%    the cantilever beam problem.
%
%    written in January 1999 by
%
%    Krister Svanberg (krille@math.kth.se)
%    Optimization and Systems Theory, KTH,
%    SE-10044 Stockholm, Sweden.
%
itte = 0;
while itte < maxite
  iter = iter+1
  itte = itte+1;

  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);

  xold2 = xold1;
  xold1 = xval;
  xval = xmma;

  [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = beam(xval);

  outvector = [f0val fval xval']'
end
