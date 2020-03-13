function [x,lambda,basis,status,it] = alps(c,A,b,l,u,sense,maxit,printlvl)
%
% alps - a linear program solver
% 
% alps assumes problem on the form 
%
%  min  z = c'x
%  s.t.  Ax | b
%        l<=x<=u
%  where | is any of =, <=, >=
%
% Calling syntax
% [x,lambda,basis,status,it] = alps(c,A,b,l,u,sense,maxit,printlvl);
% 
% input arguments:  c  cost vector
%                   A  coefficient matrix
%                   b  right-hand side
%                   l  lower bounds on x (default: 0)
%                   u  upper bound on x  (default: inf)
%               sense  sense(i)=0 means | ->  = in row i (default)
%                      sense(i)=1 means | -> <=   -"-
%                      sense(i)=2 means | -> >=   -"-
%               maxit  iteration limit (default: 10000)
%            printlvl  if > 1, iteration information is printed
%                      
% output arguments: x  solution if status = 0
%              lambda  multipliers
%               basis  indices of final basis (may include slacks)
%              status  = 0 optimal solution found
%                      = 1 problem infeasible
%                      = 2 problem unbounded
%                      = 3 iteration limit reached
%               it(1)  number of iterations in Phase I
%               it(2)  number of iterations in Phase II
%
%  NB: 'simple' is called as subroutine. 
%      Handles sparse matrices properly.
%      Use column vectors.
%
%  CAUTION: No anti-cycling strategy employed!
%
%  Author: Stefan Feltenmark <stefanf@math.kth.se>
%          Optimization & Systems Theory
%          Dept. of Mathematics, KTH
%          SE-100 44 Stockholm
%
%  Last modified 990125 by SF
%

% Tolerance to define e.g. zero
tol = sqrt(eps);
MAXIT = 10000;

% problem dimensions
[m,n] = size(A);

% Defaults
if nargin < 4 | isempty(l),
   l = zeros(n,1);
end
if nargin < 5 | isempty(u),
   u = inf*ones(n,1);
end
if nargin < 6 | isempty(sense),
   sense = zeros(m,1);
end
if nargin < 7 | isempty(maxit),
   maxit = MAXIT;
end
if nargin < 8 | isempty(printlvl),
   printlvl = 0;
end

% Check sizes
if any( size(c) ~= [n 1]) | ...
   any( size(b) ~= [m 1]) | ...
   any( size(l) ~= [n 1]) | ...
   any( size(u) ~= [n 1]) | ...
   any( size(sense) ~= [m 1]),
      whos c A b l u sense
      error('Wrong input data dimensions: check above');
end

% Check bounds 
if any( u-l <= -tol ),
   error('There are bounds with l > u');
end

% Convert bounds to general constraints
% Reverse negative unbounded variables
m1 = m;
n1 = n;
rev = find(l <= -inf);
revub = find( u(rev) <= inf );
reversed = [];
if ~isempty(revub),
   reversed = rev(revub);
   l(reversed) = -u(reversed);
   u(reversed) = inf*ones(size(reversed));
   c(reversed) = -c(reversed);
   A(:,reversed) =  -A(:,reversed);
end

% Replace free variables with diff. of two positive..
nub = find( l <= -inf );
nubfree = find( u(nub) >= inf );
free = nub(nubfree);
if ~isempty(free),
  % add columns for negative variables
  A = [A  -A(:,free)];    
  c = [c' -c(free)']';
  l(free) = zeros(size(free));
  l = [l ; zeros(size(free))];
  u = [u ; inf*ones(size(free))];
  n1 = n1 + length(free);
end

% Translate by l > -inf 
b = b - A*l;
u = u - l;

% now all variables have lower bound 0 
% Add upper bounds to constraint matrix
fub = find( u < inf);
lfub = length(fub);
I = speye(n1);
A1 = zeros(lfub, n1);
A1(:,fub) = I(fub,fub);
A = [A ; A1];
sense = [sense' ones(1,lfub)];
b = [b' u(fub)']';
m1 = m1 + lfub;
u(fub) = inf*ones(lfub,1);

% Add slack
eq = find( sense == 0 );
le = find( sense == 1 );
ge = find( sense == 2 );
I = speye(m1);
A = [A I(:,le) -I(:,ge)];
lle = length(le);
lge = length(ge);
c = [c ; zeros(lle,1) ; zeros(lge,1)];

% Reverse sign to get b >= 0
negind = find(b < -tol);
A(negind,:) = -A(negind,:);
b(negind) = -b(negind);

%%% now problem is on standard form
[x,z,lambda,basis,status,it] = simple(c,A,b,[],printlvl,maxit);

%%%% transform problem back to original form

% Reverse sign back
A(negind,:) = -A(negind,:);
b(negind) = -b(negind);

% strip slack variables
A = A(:,1:n1);
c = c(1:n1);
x = x(1:n1);
% strip upper bound rows
A = A(1:m,:);
b = b(1:m);
lambda = lambda(1:m);

% Add on lower bounds
x = x + l;

% Reconvert free variables
lfree = length(free);
if lfree > 0,
    x(free) = x(free) - x((n+1):(n+lfree));
end
x = x(1:n);

% Unreverse variables
x(reversed) = -x(reversed);




