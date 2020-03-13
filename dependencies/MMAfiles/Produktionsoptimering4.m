% Vi tar nu bort kravet startvärdet W0 och låter detta vara fritt. Samtidigt 
% kräver vi att WT=WO. Genom att lägga till kravet på rad 72 uppnås detta.
% lösningen kommer att bli fri från randvärlden och vi kan se att en punkt i 
% lösningen flyttar sig
clf
clear all

% Införande av diverse variabler och konstanter
Et=[20 30 25 10 12 13 20 35];
Wlast=[300 275 310 370 440 420 310 280 ];
Cp=[30 30 29];
Cq=[61 58 60];
C0=[110 90 120];
alfa=0.95;
Wl=10 ;
Wu=200;
W0=50;
WT=100;
Q=400; % Q Modelleras ej med, ty Wu-Wl<Q, alltså garanteras detta av Wl och Wu

% Punkter som definierar tillåtna arbetsområdet
P=[0 50;95 40;27 8;0 10;0 265;355 230;102 47;0 63;0 190;270 155;78 32;0 42;];


% Bestämmande av k och m, i ekvationen y=kx+m
for i=1:3
   for j=1:3
      K(i,j)=(P(4*i-3+j,2)-(P(4*i-4+j,2)))/(P(4*i-3+j,1)-(P(4*i-4+j,1)));
      M(i,j)=P(4*i-3+j,2)-P(4*i-3+j,1)*K(i,j);
   end
end
m=M(:);


% Genererande av c (1-24 beror av Pit, 25-48 beror på Qit, 49-56 konstant)
c=[];
for i=1:3
   c=[c -Et+Cp(i)*ones(1,8)];
end

for i=1:3
   c=[c Cq(i)*ones(1,8)];
end

c=[c zeros(1,8)]';


% Generering av b, villkorsvektorn, 1-72 relaterar till ekvationer för arbetsområdet,
% 73-80 definierar produktionbehovet, 81-88 definierar maximal ackumulering, 89-96 definierar 
% minimal ackumulering, 97 bestämmer slutvärdet
b=[];
for i=1:9
   b=[b m(i)*ones(1,8)];
end

b=[b Wlast Wu*ones(1,8) Wl*ones(1,8)]';


% Generering av sensevektorn, som betämmer beroendet mellan A och b
sense=[ones(1,24) 2*ones(1,24) 2*ones(1,24) zeros(1,8) ones(1,8) 2*ones(1,8)]';


% Genering av A vektorn som relaterar till b
A1=eye(8);
A2=zeros(8);
A=[];

for i=1:3
   A=[A;A1 A2 A2 -K(1,i)*A1 A2 A2 A2];
   A=[A;A2 A1 A2 A2 -K(2,i)*A1 A2 A2];
   A=[A;A2 A2 A1 A2 A2 -K(3,i)*A1 A2];
end

A=[A;A2 A2 A2 A1 A1 A1 diag(alfa*ones(1,7),-1)-A1+diag(alfa,7)];
A=[A;A2 A2 A2 A2 A2 A2 A1];
A=[A;A2 A2 A2 A2 A2 A2 A1];
A=sparse(A);

% Beräkning av x-vektorn, lösningsvektorn
[x,lamda]=alps(c,A,b,0*ones(56,1),inf*ones(56,1),sense);

% Plottar elproduktion(x-axel) mot elpris(y-axel) multiplicerat med 10
subplot(3,3,1)
for i=1:8 
	plot((x(i)+x(i+8)+x(i+16)),Et(i)*10,'o'), hold on
end

% Plottar värmeproduktion(x-axel) mot tankinnehåll(y-axel)
subplot(3,3,2)
for i=1:8 
   plot((x(i+24)+x(i+32)+x(i+40)),x(i+48),'X'), hold on
end

% Plottar värmeproduktion(x-axel) mot värmebehov(y-axel)
subplot(3,3,3)
for i=1:8 
   plot((x(i+24)+x(i+32)+x(i+40)),Wlast(i),'+'), hold on
end

% Plottar tillåtna området samt lösningarna för respektive verk
for i=1:3
   subplot(3,3,i+3)
   plot(P(4*i-3:i*4,1),P(4*i-3:i*4,2)),hold on
   plot(x(17+i*8:24+i*8),x(-7+i*8:i*8),'*')
end

subplot(3,3,7)
spy(A), hold on

% För att räkna ut förändringen i kostnaden, bestämmer vi en för-
% änringsvektor deltab. Kostnadsföränringen blir sedan lamda'*deltab
deltab=zeros(96,1);
deltab(75)=1;
forändringskostnaden=lamda'*deltab