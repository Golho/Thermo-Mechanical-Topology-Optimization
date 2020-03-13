function [x,z,u,basic,status,it] = simple(c,A,b,basic,printlvl,maxit)
%
%  SIMPLE - linear programming by the primal simplex method.
%    
%  Assumes problem on standard form:
%
%      minimize  c'x
%        s.t.   Ax = b   (b >= 0)
%               x >= 0
%
%  Calling syntax:
%  [x,z,u,basis,status,it] = simple(c,A,b,basic,printlvl,maxit)
%
%  Input:  c - cost vector
%          A - constraint matrix        
%          b - right-hand side
%          basic - variable index for initial basis (optional)
%          printlvl - if > 0, output is printed (optional)
%          maxit - iteration limit (optional)
%         
%  Output: x - primal solution
%          z - objective value
%          u - dual solution
%          basis - indices corr. to final basic solution
%
%          status = 0    optimal solution found
%                 = 1    no feasible solution found
%                 = 2    unbounded problem
%                 = 3    iteration limit reached
%
%          it(1)  - number of iterations in Phase I
%          it(2)  - number of iterations in Phase II
%
%  CAUTION: No anti-cycling strategy employed!
%  Please see documentation and comments in code for details.
%
%  Author: Stefan Feltenmark <stefanf@math.kth.se>
%          Optimization & Systems Theory
%          Dept. of Mathematics, KTH
%          SE-100 44 Stockholm
%
%  Last modified 1999/01/29 by SF
%

%
% Constants
%
tol = sqrt(eps);
MAXIT = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m1,n]=size(A);
m = m1;
if( nargin < 6)
 maxit = MAXIT;
end
%
% Check input data
%
if m > n,
  error('Overdetermined system (unless rank deficient): Use x = A\b instead!');
end %if

if length(c) ~= n | length(b) ~= m,
    error('Data of incompatible size');
end % if

if any(b < -tol),
 error('Right-hand side must be non-negative.');
end

% Default printlvl = 0 
if nargin <= 4,
 printlvl = 0;
end
%
% Check suggested start basis
%
if nargin > 3
  B = A(:,basic);
  % a bit coarse but short...
  if rank(full(B)) == m,
   Bi = inv(B);
   y0 = Bi*b;
   if y0 >= 0,
    phase = 2;
    q = zeros(n,1);
    q(basic) = ones(m,1);
    nonbasic = find(~q);
   else
    if printlvl > 0,	
     disp('Initial basis not feasible!');
    end % if
   end % if
  else
     if printlvl > 0,
      disp('Initial basis not invertible...doing phase I');
     end % if
     phase = 1;
  end % if
else
 phase = 1;
end % if

if phase == 1,
    cc = c;
    c = [zeros(n,1) ; ones(m,1)];
    A = [ A speye(m)];
    B = speye(m);
    Bi = B;
    nonbasic = (1:n)';
    basic = n+1:n+m';
    y0 = b;
    z = c(basic)'*y0;
end % if

rows = 1:m; % number of rows may decrease

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%  Iterate until some terminiation criteria satisfied   %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it  = 0;
it1 = 0;
it2 = 0;
while 1,
%
% Same loop both phases
%
 if printlvl > 0,
  disp(sprintf('Phase %d',phase));
  disp(sprintf(' It   out   in    red. cost         step    objective'));
 end %
 
 rowcnt = 0;
 while it < maxit,
    
 it = it + 1;
 if printlvl > 0 & rowcnt == 20,
  disp(sprintf(' It   out   in    red. cost         step    objective'));
  rowcnt = 0;
 end %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   Compute dual solution and reduced costs             %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  u = c(basic)'*Bi;
  r = c(nonbasic)' - u*A(rows,nonbasic);
  
  % Leaving variable index
  i = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   Determine entering variable                         %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nind = find(r < -tol);
  if ~isempty(nind),

    % Dantzig's rule
    [rj,j]=min(r);
    j = j(1);

    % Steepest edge pricing
    % Y = Bi*A(rows,nonbasic);
    % nd = sqrt(1+diag( Y'*Y ));	
    % [rj,j]=min(r./nd);
    % j = j(1);

  elseif phase == 1,

    % Infeasibility left? Then (primal) problem infeasible.
    if( z > tol )
      flag = 1;
      break;
    end 

    % Are there artificial variables left?
    artind = find( basic > n );
    if ~isempty( artind ),
%
% Find original variables that can be pivoted into the basis
% or conclude that there is redundant equations.
%
     ind = find( nonbasic' <= n );
     if isempty( ind ),
        % This cannot happen since we check m <= n 
        flag = 1;
        break;
      else      
        %
        % Check if we can replace artificial variables
        % with original ones, and still maintain a bfs.
        %
        Y = Bi*A(rows,nonbasic(ind));
        Y = Y(artind,:);
        j = 0;
        ndel = 0;
        for k = 1:length( artind ),
         for l = 1:length(ind),
           if abs(Y(k,l)) > tol,
              % Choose l'th nonbasic to enter
              j = ind(l);
              % Choose k'th basic artificial to leave
              i = artind(k);
              break;
           end % if
         end % for l
         if j > 0,
            break;
         else
           % We can cross out a row!
           col = artind(k);
           row = basic(col) - n;
           % actually col == row!!

           if printlvl > 0,
            disp(sprintf('Linear dependence detected: deleting row %d',row));
           end

           % Delete row index	
           rows(row-ndel:m-1) = rows(row-ndel+1:m);
           rows = rows(1:m-1);

           % Delete from list of basics	
           basic(col:m-1) = basic(col+1:m);
           basic = basic(1:m-1);

           % We reinvert for now. Simple but not very efficient.
           B = A(rows, basic);
           Bi = inv(B);
           
           % position of artificial variables in basic is shifted
           artind = artind - 1;

           % Delete from y0
           y0(col:m-1) = y0(col+1:m);
           y0 = y0(1:m-1);

           m = m - 1;
           ndel = ndel + 1;
         end % if
        end % for k
        if j == 0,
          flag = 0;
          break;
        end % if
       end % if ind
    else
      % no artificial variables left in basis!
      flag = 0;
      break;
    end % if artind
  else
    % optimal solution found!
    flag = 0;
    break;
  end % if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%     Determine leaving variable                        %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  yj = Bi*A(rows,nonbasic(j));
  if i == 0,

    pind = find(yj > tol);
    if isempty( pind ),
    % Problem unbounded!
     flag = 2;
     break;
    end % if
%  
% Max steplength, and basics that hit zero
%
    [alpha,i] = min( y0(pind)./yj(pind) );
    i = pind(i);
    i = i(1);

  else
   alpha = 0;
  end % if
  
  % Update solution.
  y0 = y0 - alpha*yj;
  y0(i) = alpha;
  
  % Update basis.
  temp = basic(i);
  basic(i) = nonbasic(j);
  nonbasic(j) = temp;
  M = speye(size(Bi));
  M(:,i) = -yj./yj(i);
  M(i,i) = 1/yj(i);
  % Pivot
  Bi = M*Bi;

  z = c(basic)'*y0;
  if printlvl > 0,
   disp(sprintf('%4.0d %4.0d %4.0d %12.3f %12.3f %12.3f', ...
                              it,nonbasic(j),basic(i),rj,alpha,z));
   rowcnt = rowcnt + 1;
  end % if


 end % while
 
  if it >= maxit,
    flag = 3;
  end % if

 if flag == 0,
  if phase == 1,
    it1 = it; 
    phase = 2;
    c = cc;
    q = zeros(n,1);
    q(basic) = ones(m,1);
    nonbasic = find(~q);
 else
    it2 = it - it1;
    if printlvl > 0,
     disp('Solution is optimal.'); 
    end % if
    x = zeros(n,1);
    x(basic) = y0;
    z = c'*x;
    status = 0;
    break;
  end % if
 elseif flag == 1,
    if printlvl > 0,
     disp('Problem infeasible.');
    end % if
    x = zeros(n,1);
    x(basic) = y0;
    z = Inf;
    status = 1;
    it1 = it;
    break;
 elseif flag == 2,
    if printlvl > 0,
     disp('Problem unbounded.');
    end % if
    x = zeros(n,1);
    x(basic) = yj;
    x(j) = 1;
    z = -Inf;
    status = 2;
    it2 = it-it1;
    break;
 elseif flag == 3,
    if printlvl > 0,
     disp('Iteration limit reached.');
    end % if
    x = zeros(n,1);
    if phase == 1,
       it1 = it;
       ind = find(basic <= n);
       x(basic(ind)) = y0(ind);
       z = c(basic)'*y0;
       basic = basic(ind);
    else
      x(basic) = y0;
      z = c'*x;
      it2 = it-it1;
    end
    status = 3;
    break;
 end % if
end % while

% assign output arguments
it = [it1 it2]';

% consider case when rows deleted
u1 = zeros(m1,1);
u1(rows) = u;
u = u1;


