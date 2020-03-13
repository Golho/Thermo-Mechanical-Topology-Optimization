function [x,ub,lb,bbtree,status] = lips( c,A,b,l,u,sense,I,opttol,maxit,bs,printlvl )
%
%  LIPS.M - Linear Integer Program Solver
%  Method: branch-and-bound with LP relaxations.                 
%  The following form of the problem is assumed:
%
%  Minimize   c'x
%
%   s.t.    Ax | b,
%           l<= x <=u,
%           x(i), i in I integer      
%
% where | is one of =, >=, or <= for each row.
%
% Calling syntax:
% [x,ub,lb,bbtree,status] = lips( c,A,b,l,u,sense,I,opttol,maxit,bs,printlvl )
%
% Input  c  cost vector
%        A  constraint matrix
%        b  right-hand side
%        l  variable lower bound (default: 0)
%        u  variable upper bound (default: 1)
%    sense  sense(i)=0 means | ->  = in row i (default)
%           sense(i)=1 means | -> <=   -"-
%           sense(i)=2 means | -> >=   -"-
%        I  index vector for integer variables (default: all)
%   opttol  optimality tolerance (default: 1.0e-6)
%    maxit  iteration limit (default: 1000)
%       bs  branching strategy: 1 = depth first
%                               2 = bredth first
%                               3 = best first (default)
% printlvl  if printlvl >= 1, one line per iteration is printed
%
%   Output
%        x  best solution found
%       ub  best known upper bound
%       lb  least lower bound
%   bbtree  matrix representing b&b-tree
%           each row corr. to node, the columns mean
%           bbtree(:,1) node number 
%           bbtree(:,2) parent node number 
%           bbtree(:,3) status: 0 = active, 1 = infeasible, 2 = fathomed,
%                               3 = dominated, 4 = relaxed soln. feasible
%           bbtree(:,4) depth of node in tree (root has depth 0)
%           bbtree(:,5) local lower bound
%           bbtree(:,6) local upper bound
%           bbtree(:,7) 0 = left (<=) child, 1 = right (>=) child.
%           bbtree(:,8) index of branching variable
%           bbtree(:,9) branching variable bound 
%   status  = 0, solution optimal within tolerance
%           = 1, feasible solution found
%           = 2, no feasible solution found
%           = 3, problem infeasible
%           
%
%  NB: 'lips' uses 'alps' to solve the LP-relaxations. 
%       Handles sparse matrices properly.
%       Use column vectors.
% 
%
%  Author: Stefan Feltenmark <stefanf@math.kth.se>
%          Optimization & Systems Theory
%          Dept. of Mathematics, KTH
%          SE-100 44 Stockholm
%                                        
%                                                     (1999/01/30)
%


%%%%%%%%%%%%%%%%%%%%%%%%%%  Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPTTOL = 1.0e-6;
FEATOL = 1.0e-5;
MAXIT  = 1000;
BSTRAT = 3; % Branching strategy: 1 = depth first
            %                     2 = bredth first
            %                     3 = best first (recommended)


%%%%%%%%%%%% No parameters to be changed below this line %%%%%%%%%%%%
[m,n] = size(A);
xstar   = NaN*ones(n,1);
ub      = inf;
lb      = -inf;
nactive = 0;
active  = [];
bbtree  = [];
bbpath  = [];
glbs = [];

% Defaults
if nargin < 4 | isempty(l),
   l = zeros(n,1);
end
if nargin < 5 | isempty(u),
   u = ones(n,1);
end
if nargin < 6 | isempty(sense),
   sense = zeros(m,1);
end
if nargin < 7 | isempty(I),
   I = (1:n)';
end
if nargin < 8 | isempty(opttol),
   opttol = OPTTOL;
end
if nargin < 9 | isempty(maxit),
   maxit = MAXIT;
end
if nargin < 10 | isempty(bs),
   bs = BSTRAT;
end
if nargin < 11 | isempty(printlvl),
   printlvl = 0;
end
%%% beginning of solution procedure
it = 1;

if printlvl > 0,
   disp(sprintf('  It #act curr        lb       glb        ub       gub       gap st dp '));
end

% Index of =, >=, ,= rows
eq = find(sense == 0);
le = find(sense == 1);
ge = find(sense == 2);

% add root node to tree
nactive         = 1;
active          = [1];
depth           = 0;
nodecnt         = 1;
bbtree          = [bbtree ; nodecnt  0  0  0 -inf inf 0 NaN NaN];

% initialize bounds et c
ub		 =  inf;
lb     = -inf;
gap	 =  inf;
done   =  0;
rowcnt =  0;
status = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  MAIN LOOP                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while it <= maxit & ~done,

 % choose an active node
 % [current] = selectnode(bs);
 %%
 active = find(bbtree(:,3) == 0);
 nact = length(active);

 if nact == 0,
    current = 0;
    done = 1;
    break;
 end

 if bs == 1,
    %% depth first, last one created
    current = active(nact);
 elseif bs == 2,
    %% bredth first, the active with lowest depth
    [ignore,mind] = min(bbtree(active,4));
    current = active(mind);
 elseif bs == 3,
    %% best first, i.e any with lowest lower bound
    [ignore,minl] = min(bbtree(active,5));
    current = active(minl);
 else
   error('Parameter bs has value out of range');
 end

 if current == 0,
   % No active node left
   done = 1;
   break;
 end

 depth =  bbtree(current,4);

 % include parent node branchings
 % trace path to row
end