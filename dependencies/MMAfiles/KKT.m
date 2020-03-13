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
round((lam.*s)*1e7)

disp('KKT-villkor 3.2k')
round((xsi.*(xmma-xmin))*1e7)

disp('KKT-villkor 3.2l')
round((eta.*(xmax-xmma))*1e7)

disp('KKT-villkor 3.2m')
round((mu.*ymma)*1e7)

disp('KKT-villkor 3.2n')
round((zet.*zmma)*1e7)
