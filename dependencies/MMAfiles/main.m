clf

maxite=20;
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

  [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = pro(xval);

  outvector = [f0val fval' xval']'
end

p = [0 0 0 0 0 0 0 -1000]';
q = [0 0 0 0 0 -1000 0 0]';
R = 760518/3200*[1 0 .5 0 0 -1 -.5 0 0 0; 
		 0 0 .5 0 1 0 .5 0 0 0;
		 0 .5 0 1 0 0 0 -.5 -1 0;
		 0 -.5 0 0 -1 0 0 -.5 0 0;
		 0 0 0 0 0 1 0 .5 0 0;
		 0 0 0 0 0 0 0 .5 0 1;
		 0 0 0 0 0 0 .5 0 1 0;
		 0 0 0 0 0 0 -.5 0 0 -1];
Kinv=inv(R*diag(xval)*R');

u=Kinv*p;uu=[];
v=Kinv*q;vv=[];

xy=[0 0 1 1 2 2; 1 0 1 0 1 0]';
andpkter=[1 3;1 4;2 3;2 4;3 4;3 5;3 6;4 5;4 6;5 6];


plot(xy(1,:),xy(2,:),'o'), axis equal, hold on

for i=1:10
l=line([xy(andpkter(i,1),1) xy(andpkter(i,2),1)],[xy(andpkter(i,1),2) xy(andpkter(i,2),2)]);
%set(l,'LineWidth',xval(i)/4);
end

for i=1:4
uu=[uu; u(i*2-1) u(i*2)];
vv=[vv; v(i*2-1) v(i*2)];
end

xy1=xy+[zeros(2,2); 10*uu];

for i=1:10
l=line([xy1(andpkter(i,1),1) xy1(andpkter(i,2),1)],[xy1(andpkter(i,1),2) xy1(andpkter(i,2),2)]);
end

xy2=xy+[zeros(2,2); 10*vv];

for i=1:10
l=line([xy2(andpkter(i,1),1) xy2(andpkter(i,2),1)],[xy2(andpkter(i,1),2) xy2(andpkter(i,2),2)]);
end

title('lastfall 1 & 2, förskjutningar förstärkta tio gånger')

%kontroll KKT-villkor


disp('KKT-villkor 3.2a')
round ( (df0dx+(lam'*dfdx)'-xsi+eta) *1e7 )
disp('KKT-villkor 3.2b')
round( ( c+d.*ymma-lam-mu ) *1e7)

disp('KKT-villkor 3.2c')
round((a0-zet-lam'*a)*1e7)
disp('KKT-villkor 3.2d')
round((fval-a.*zmma-ymma+s)*1e7)

disp('KKT-villkor 3.2e')
xmma
disp('KKT-villkor 3.2f')
zmma  
ymma  
s     

disp('KKT-villkor 3.2g')
xsi
eta

disp('KKT-villkor 3.2h')
zet
mu

disp('KKT-villkor 3.2i')
lam

disp('KKT-villkor 3.2j')
lam.*s

disp('KKT-villkor 3.2k')
xsi.*(xmma-xmin)

disp('KKT-villkor 3.2l')
eta.*(xmax-xmma)

disp('KKT-villkor 3.2m')
mu.*ymma

disp('KKT-villkor 3.2n')
zet.*zmma
