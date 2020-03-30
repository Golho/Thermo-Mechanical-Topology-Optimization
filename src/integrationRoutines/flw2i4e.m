function [Ke]=flw2i4e(ex,ey,ep,D)
% Ke=flw2i4e(ex,ey,ep,D)
% [Ke,fe]=flw2i4e(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness (conductivity)
%  matrix for 4 node isoparametric field element
%
% INPUT:  ex = [x1 x2 x3 x4]   element coordinates
%         ey = [y1 y2 y3 y4]
%                             
%         ep = [t ir]          thickness and 
%                              integration rule
%
%         D  = [kxx kxy;
%               kyx kyy]       constitutive matrix
%
%         eq                   heat supply per unit volume
%
% OUTPUT: Ke :  element 'stiffness' matrix (4 x 4)
%         fe :  element load vector (4 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-24
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  t=ep(1); ir=ep(2); ngp=ir*ir;

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


  Ke1 = zeros(4,4);
  B = zeros(2, 4);
  JT=dNr*[ex;ey]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=abs(det(JT(indx,:)));
    B = JT(indx, :) \ dNr(indx,:);
    Ke1=Ke1+B'*D*B*detJ*wp(i);
  end

  Ke=Ke1*t;
%--------------------------end--------------------------------
