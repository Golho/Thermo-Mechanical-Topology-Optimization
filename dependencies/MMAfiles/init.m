  m = 2;
  n = 10;
  xval  = 15*ones(n,1);
  xold1 = zeros(n,1);
  xold2 = zeros(n,1);
  low   = zeros(n,1);
  upp   = zeros(n,1);
  xmin  = ones(n,1);
  xmax  = 20*ones(n,1);
  c = 1000*ones(m,1);
  d = zeros(m,1);
  a0 = 1;
  a = zeros(m,1);
  [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = pro(xval);
  
  outvector = [f0val fval' xval']'
  iter = 0;
