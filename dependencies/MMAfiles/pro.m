function [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = pro(x);


f0val = ones(1,10)*x
df0dx = ones(10,1);
df0dx2 = zeros(10,1);

X=diag(x);
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

K = R*X*R';
Kinv = inv(K);

fval = [p'*Kinv*p-10; q'*Kinv*q-10];

dfdx=[];
dfdx2=[];

for i=1:10
  A1=Kinv*R(:,i)*R(:,i)'*Kinv;
  A2=A1*R(:,i)*R(:,i)'*Kinv;
  dfdx=[dfdx -[p'*A1*p; q'*A1*q]];
  dfdx2=[dfdx2 2*[p'*A2*p; q'*A2*q]];
end



