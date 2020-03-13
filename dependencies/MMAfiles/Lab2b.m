clf
clear all

% Införande av diverse variabler och konstanter
K=[50 50 55 75 75 70 65 60];
C=3000;
rabatt=20; %Anges i procent 
D=C*(100-rabatt)/100;
L=100;
r=100;
s=800; % Godtyckligt stort tal
Zmax=80;

% Genererande av c 
c=[L*ones(1,8) C*ones(1,8) D*ones(1,8) zeros(1,8)]';


% Generering av b
b=[K zeros(1,8) zeros(1,8)]';

% Generering av sensevektorn, som betämmer beroendet mellan A och b
sense=[zeros(1,8) 2*ones(1,8) ones(1,8) ]';


% Genering av A vektorn som relaterar till b
A1=eye(8);
A2=zeros(8);

A=[diag(ones(1,7),-1)-A1 A1 A1 A2];
A=[A;A2 A1 A2 -r*A1];
A=[A;A2 A2 A1 -(s-r)*A1];
A=sparse(A);

% Genering av u & l & I

I=(25:32);
u=[80*ones(1,8) r*ones(1,8) (s-r)*ones(1,8) ones(1,8)]';
l=[zeros(1,32)]';

[x]=lips(c,A,b,l,u,sense,I);
x(9:16)
x(17:24)

% Utskrift av optimalt inköpsprogram
inkop=[];
vecka=[];
for i=1:8
  inkop = [inkop x(8+i)+x(16+i)];
  vecka = [vecka i];
end
plot(vecka,inkop,'mx-'), hold on

% Beräkning av minskad kostnad m.m.
optimeradkostnad=c'*x;
fortjanst=1500000-optimeradkostnad
konsultarvode = 2 * 8 * 1000