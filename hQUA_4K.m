function [Ke] = hQUA_4K(ex, ey, D, t)
%FLW2RE Compute the element heat stiffness matrix for a 2D Melosh element
% It is assumed that the conductivity matrix is symmetric.
%   ex      x-coordinates of the nodes [m]
%   ey      y-coordinates of the nodes [m]
%   D       Thermal conductivity matrix [2 x 2]
%   t       Thickness [m]

t=ep(1); ir=ep(2); ngp=ir*ir;
  if nargin==4; eq=0 ; end

  if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
  elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
  elseif ir==3
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
  else
    disp('Used number of integration points not implemented');
    return
  end
  wp=w(:,1).*w(:,2);


  xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;


  dNr(1:2:r2,1)=-(1-eta)/4;     dNr(1:2:r2,2)= (1-eta)/4;
  dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
  dNr(2:2:r2+1,1)=-(1-xsi)/4;   dNr(2:2:r2+1,2)=-(1+xsi)/4;
  dNr(2:2:r2+1,3)= (1+xsi)/4;   dNr(2:2:r2+1,4)= (1-xsi)/4;


  Ke1=zeros(4,4);  fe1=zeros(4,1);
  JT=dNr*[ex;ey]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminanten lika med noll!')
    end
    JTinv=inv(JT(indx,:));
    B=JTinv*dNr(indx,:);
    Ke1=Ke1+B'*D*B*detJ*wp(i);
  end

  Ke=Ke1*t;  fe=fe1*t*eq;




area = abs(ex(3)-ex(1))*abs(ey(3)-ey(1));

s1 = -D(1, 1)/6  -    D(2, 2)/6;
s2 = -D(1, 1)/3 +    D(2, 2)/6;
s3 = D(1, 1)/6 -    D(2, 2)/3;
s4 = D(1, 1)/3 +    D(2, 2)/3;

I = [s4,     s2,     s1,     s3;
     s2,     s4      s3      s1;
     s1      s3      s4      s2;
     s3      s1      s2      s4];
I2 = [4,     -1,     -2,     -1;
      -1,     4      -1      -2;
      -2      -1      4      -1;
      -1      -2      -1      4] * 1/6;
Ke = t * area*I2;
end
