%    This is the file beaminit.m
%    in which some vectors for the cantilever
%    beam problem are initialized.
%
%    written in January 1999 by
%
%    Krister Svanberg (krille@math.kth.se)
%    Optimization and Systems Theory, KTH,
%    SE-10044 Stockholm, Sweden.
%
  m = 1;
  n = 5;
  xval  = 5*ones(n,1);
  xold1 = zeros(n,1);
  xold2 = zeros(n,1);
  low   = zeros(n,1);
  upp   = zeros(n,1);
  xmin  = ones(n,1);
  xmax  = 10*ones(n,1);
  c = 1000*ones(m,1);
  d = zeros(m,1);
  a0 = 1;
  a = zeros(m,1);
  [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = beam(xval);
  outvector = [f0val fval xval']'
  iter = 0;
