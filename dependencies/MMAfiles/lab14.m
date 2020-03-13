clf
clear all

et=[20 30 25 10 12 13 20 35];
Wlast=[300 275 310 370 440 420 310 280];
c0=[110 90 120];
cp=[30 30 29];
cq=[61 58 60];
alfa=.95;
Wl=10;
Wu=200;
W0=50;
WT=100;
Q=400;


% tillåtna området

P=[0 50;95 40;27 8;0 10;0 265;355 230;102 47;0 63;0 190;270 155;78 32;0 42];



% ekvationerna y=kx+m för begränsningarna

for i=1:3

  for j=1:3

    K(i,j)=(P(i*4-4+j,2)-P(i*4-3+j,2))/(P(i*4-4+j,1)-P(i*4-3+j,1));

    M(i,j)=P(i*4-4+j,2)-K(i,j)*P(i*4-4+j,1);

  end

end

m=M(:);

% generering av b, vilkorsvektor

b=[];

for i=1:9
b=[b ones(1,8)*m(i)];
end

b=[b Wlast Wu*ones(1,8) Wl*ones(1,8)];
b=b';

% generering av c
c=[];
for i=1:3
c=[c cp(i)*ones(1,8)-et];
end
for i=1:3
c=[c cq(i)*ones(1,8)];
end
c=[c zeros(1,8)]';

% generering av sense

sense=[ones(1,24) 2*ones(1,24) 2*ones(1,24) zeros(1,8) ones(1,8) 2*ones(1,8)]'; 

% generering av A

A1=eye(8);
A2=zeros(8);

A=[];

for i=1:3
A=[A;A1 A2 A2 -K(1,i)*A1 A2 A2 A2];
A=[A;A2 A1 A2 A2 -K(2,i)*A1 A2 A2];
A=[A;A2 A2 A1 A2 A2 -K(3,i)*A1 A2];
end


A=[A;A2 A2 A2 A1 A1 A1 diag(alfa*ones(1,7),-1)-A1+diag(alfa,7)];
A=[A;zeros(8,48) A1];
A=[A;zeros(8,48) A1];
A=sparse(A);

x=alps(c,A,b,0*ones(56,1),inf*ones(56,1),sense);


% rita tillåtna området och lösningarna

for i=1:3

subplot(2,2,i)
plot(P(i*4-3:i*4,1),P(i*4-3:i*4,2)), hold on

plot(x(17+i*8:24+i*8),x(-7+i*8:i*8),'o'), hold on
plot(x(19+i*8),x(-5+i*8),'x'), hold on
end

subplot(2,2,4)
spy(A)

[A*x b A*x-b<=0]

qpW=[];

for i=1:3
qpW=[qpW x(i*8+17:i*8+24) x(i*8-7:i*8)];
end
qpW=[qpW x(49:56)]

branslekostnad=c'*x+c0*ones(3,1)*8




