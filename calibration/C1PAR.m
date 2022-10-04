
function [Y1c, Tta] = C1PAR(x, Prss, Pctt)

global frequ;
global kp;

Bcnt=x(1);
rhob=x(2).*10.^kp;

mm = 1e-3;
pico =1e-12;

LC=Pctt(3).*mm;

lmbd = Pctt(6);

ms = Pctt(2)*1e-3;

% propriedades cristalograficas do LaNi5
dLaNi5 = 8320.60623983755; %kg/m3

% volume do solido pela densidade (m3)
Vs = (ms./dLaNi5); % ms é a massa da amostra

% volume do anel total ( Lc aLcura total da amostra)
Vt = pi.*LC.*(((7.5001.*mm).^2) - ((5.7.*mm).^2)); %m3

% porosidade
Tta = (Vt-Vs)./Vt;

omega=2.*pi.*frequ;

LAMB=lmbd.*mm;

% constantes dieletricas
epsilon0 = 8.854187817*pico; % vacuo
% constante dieletrica é " e=er*e0 "
EQ = 4.2.*epsilon0; % permitividade relativa quartzo
eAe0=1.00007.*epsilon0;%1.00007;
eBe0=2.4211E-4.*epsilon0;% 2,45495E-4
% resistividade do cobre
rhocu= 1.68*(10^-8); 

% geometria fixa
Ld = 7.09.*mm; % altura da base do eletrodo interno
SGauss = 37.*mm; %altura total gaussiana
Gv = SGauss-Ld; % altura variavel da gaussiana

r1 = 5.7.*mm;
r2 = 7.4999.*mm;
r3 = 7.5001.*mm;
r4 = 9.0*mm;
r5 = 9.15*mm;
r6 = 5.925*mm;
r7 = 7.475*mm;

Cut = (log(r6./r1))./(2.*pi);
Cu = (log(r7./r6))./(2.*pi);

% constantes geometrica simetria axial
A = (2.*pi)./(log(r3./r7));
Dd = (2.*pi)./(log(r3./r2));
E = (2.*pi.*EQ)./(log(r4./r3));
F = (2.*pi)./(log(r5./r4));

alpha1 = Cut + Cu;
alpha2 = E.*F;

Hs = (log(r3./r1))./(2.*pi);

Hp = (2.*pi)./(log(r3./r1));

% Correções da camara

  Ap1= 1.5909e-07;
  Ap2= 1.0275e-05;
  Ap3= -0.10297;
  Bp1= -7.0222e-07;
  Bp2= -3.9514e-05;
  Bp3= 1.1678;
  Cp1= -1.2967e-06;
  Cp2= 0.00012059;
  Cp3= -2.369;
  
Apc = Ap1.*Prss.^2+Ap2.*Prss+Ap3;
Bpc = Bp1.*Prss.^2+Bp2.*Prss+Bp3;
Cpc = Cp1.*Prss.^2+Cp2.*Prss+Cp3;

Zcf1jx = (Apc.*lmbd.^2+Bpc.*lmbd+Cpc).*pico;

  ar = 681443.7497;
  br = 5801399999999.896;
  cr = 0.22611;

% Zrf1ix = ar + br.*exp(-lmbd./cr);
% ponto de maximo em lambda 4mm 1.568e6
% ponto de minimo em lambda 4mm 5.685e5
% deta de +/- 0.54975e6
Zrf1jx = (ar + br.*exp(-lmbd./cr))-0.54975e6;

rcnt = 1;

pcnt = Bcnt;

% função total modelo Zlambda
ZTTL=(1j).*(-A.*(Cut.^2).*rhocu.^2.*alpha2.*((eAe0+eBe0.*Prss).^2).*omega.^2+(1j).*Cut.*rhocu.*((1+(-2.*Cut+2.*alpha1).*A).*alpha2+A.*((eAe0+eBe0.*Prss).*F+E)).*(eAe0+eBe0.*Prss).*omega-((alpha2+A.*((eAe0+eBe0.*Prss).*F+E)).*(Cut-alpha1))).*(rcnt.*(-LAMB+LC).*rhob.*pcnt.*Hs.*omega.*(eAe0+eBe0.*Prss).*((eAe0+eBe0.*Prss).*F+E).*Hp.^2+(-E.*F.*Hs.*LAMB.*rcnt.*rhob.*pcnt.*(eAe0+eBe0.*Prss).*omega+(1j).*(LAMB-LC+LC.*pcnt).*((eAe0+eBe0.*Prss).*F+E)).*Hp+(1j).*E.*LAMB.*F).*(-Zrf1jx.*omega.*Zcf1jx+(1j)).*(alpha2+(eAe0+eBe0.*Prss).*Dd.*F+Dd.*E)./omega./((-(Cut.^2).*Zrf1jx.*Zcf1jx.*A.*alpha2.*rhocu.^2.*((eAe0+eBe0.*Prss).^3).*((((Ld.*(eAe0+eBe0.*Prss).*Dd)+E.*LAMB).*F+(E.*Dd.*Ld)).*alpha2+E.*LAMB.*F.*Dd.*((eAe0+eBe0.*Prss).*F+E)).*omega.^3+(1j).*Cut.*rhocu.*(((((Ld.*(eAe0+eBe0.*Prss).*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*Dd+E.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*LAMB+Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*Cut+2.*((Ld.*alpha1.*(eAe0+eBe0.*Prss).*Dd)+(-eAe0./0.2e1+(E.*alpha1)-(eBe0.*Prss)./0.2e1).*LAMB+((Gv.*(eAe0+eBe0.*Prss))./0.2e1)).*Zrf1jx.*Zcf1jx).*F+((Ld.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*Dd+Zcf1jx.*rhocu).*Cut+2.*Zrf1jx.*Zcf1jx.*((Dd.*Ld.*alpha1)-LAMB./2+(Gv./0.2e1))).*E).*A+Zrf1jx.*Zcf1jx.*(((Ld.*(eAe0+eBe0.*Prss).*Dd)+E.*LAMB).*F+(E.*Dd.*Ld))).*(alpha2.^2)+(((Dd.*(E.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*LAMB+Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*Cut+2.*Zrf1jx.*(((-eAe0./0.2e1+(E.*alpha1)-(eBe0.*Prss)./0.2e1).*LAMB+(((eAe0+eBe0.*Prss).*(Ld+Gv))./0.2e1)).*Dd+E.*LAMB./2).*Zcf1jx).*F+(Cut.*rhocu+Zrf1jx.*(Gv-LAMB+Ld)).*Zcf1jx.*E.*Dd).*A+E.*F.*Dd.*Zcf1jx.*Zrf1jx.*LAMB).*((eAe0+eBe0.*Prss).*F+E).*alpha2+Zrf1jx.*Zcf1jx.*E.*A.*LAMB.*F.*Dd.*(((eAe0+eBe0.*Prss).*F+E).^2)).*((eAe0+eBe0.*Prss).^2).*omega.^2+(((-2.*(rhocu.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx+E.*LAMB).*(Cut.^2)+(-Ld.*alpha1.*rhocu.*(eAe0+eBe0.*Prss).*Dd+(rhocu.*eBe0.*Prss./2+rhocu.*eAe0./2-E.*rhocu.*alpha1-Zrf1jx.*Zcf1jx./2).*LAMB-Gv.*Prss.*eBe0.*rhocu./2-Gv.*eAe0.*rhocu./2+Zcf1jx.*(-2.*alpha1.*rhocu+Gv.*Zrf1jx)./2).*Cut-Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB)./2).*(eAe0+eBe0.*Prss).*F-2.*E.*(rhocu.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx).*(Cut.^2)+(-Ld.*alpha1.*rhocu.*(eAe0+eBe0.*Prss).*Dd+(rhocu.*eBe0.*Prss./2+rhocu.*eAe0./2-Zrf1jx.*Zcf1jx./2).*LAMB-Gv.*Prss.*eBe0.*rhocu./2-Gv.*eAe0.*rhocu./2+Zcf1jx.*(-2.*alpha1.*rhocu+Gv.*Zrf1jx)./2).*Cut-Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB)./2)).*A+((Ld.*(eAe0+eBe0.*Prss).*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*Dd+E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB+Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*Cut+((Ld.*(eAe0+eBe0.*Prss).*Dd)+E.*LAMB).*Zrf1jx.*Zcf1jx.*alpha1).*F+((Ld.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*Dd+Zcf1jx.*rhocu).*Cut+Dd.*Ld.*Zcf1jx.*Zrf1jx.*alpha1).*E).*(alpha2.^2)-2.*(((Dd.*rhocu.*(eAe0+eBe0.*Prss).*(Zcf1jx+E.*LAMB).*(Cut.^2)+(-((-rhocu.*eBe0.*Prss./2-rhocu.*eAe0./2+E.*rhocu.*alpha1+Zrf1jx.*Zcf1jx./2).*LAMB+rhocu.*eBe0.*(Ld+Gv).*Prss./2+rhocu.*(Ld+Gv).*eAe0./2-Zcf1jx.*(-2.*alpha1.*rhocu+Zrf1jx.*(Ld+Gv))./2).*(eAe0+eBe0.*Prss).*Dd-E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB./2-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)./2).*Cut-((eAe0+eBe0.*Prss).*(Gv-LAMB+Ld).*Dd+E.*LAMB).*Zrf1jx.*Zcf1jx.*alpha1./2).*F-(-2.*rhocu.*Zcf1jx.*Dd.*(Cut.^2)+(((-rhocu.*eBe0.*Prss-rhocu.*eAe0+Zrf1jx.*Zcf1jx).*LAMB+rhocu.*eBe0.*(Ld+Gv).*Prss+rhocu.*(Ld+Gv).*eAe0-Zcf1jx.*(-2.*alpha1.*rhocu+Zrf1jx.*(Ld+Gv))).*Dd+Zcf1jx.*rhocu).*Cut+Dd.*Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB+Ld)).*E./2).*A-(((E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB+Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*Cut+E.*Zcf1jx.*Zrf1jx.*alpha1.*LAMB).*F+E.*Cut.*Zcf1jx.*rhocu).*Dd./2).*((eAe0+eBe0.*Prss).*F+E).*alpha2+(((E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB+Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*Cut+E.*Zcf1jx.*Zrf1jx.*alpha1.*LAMB).*F+E.*Cut.*Zcf1jx.*rhocu).*A.*Dd.*(((eAe0+eBe0.*Prss).*F+E).^2)).*(eAe0+eBe0.*Prss).*omega+(1j).*(Cut-alpha1).*(((Gv-LAMB).*(eAe0+eBe0.*Prss).*((eAe0+eBe0.*Prss).*F+E).*A+((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx+E.*LAMB).*(eAe0+eBe0.*Prss).*F+E.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx)).*(alpha2.^2)+((((eAe0+eBe0.*Prss).*(Gv-LAMB+Ld).*Dd+Zcf1jx+E.*LAMB).*(eAe0+eBe0.*Prss).*F+((eAe0+eBe0.*Prss).*(Gv-LAMB+Ld).*Dd+Zcf1jx).*E).*A+((eAe0+eBe0.*Prss).*(Zcf1jx+E.*LAMB).*F+E.*Zcf1jx).*Dd).*((eAe0+eBe0.*Prss).*F+E).*alpha2+((eAe0+eBe0.*Prss).*(Zcf1jx+E.*LAMB).*F+E.*Zcf1jx).*A.*Dd.*(((eAe0+eBe0.*Prss).*F+E).^2))).*rcnt.*(-LAMB+LC).*rhob.*pcnt.*Hs.*omega.*(eAe0+eBe0.*Prss).*Hp.^2+(A.*E.*F.*Dd.*Hs.*Ld.*(Cut.^2).*pcnt.*rcnt.*rhob.*Zcf1jx.*Zrf1jx.*rhocu.^2.*(alpha2.^2).*LAMB.*((eAe0+eBe0.*Prss).^4).*omega.^4+(-1.*1j).*Cut.*alpha2.*rhocu.*(((((((Hs.*pcnt.*rcnt.*rhob.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*E+(eAe0+eBe0.*Prss).*Zrf1jx.*Zcf1jx.*rhocu).*LAMB+LC.*Zrf1jx.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Ld.*Dd+(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB.*rhocu).*Cut+2.*rcnt.*rhob.*Zrf1jx.*pcnt.*Zcf1jx.*E.*((Dd.*Ld.*alpha1)-LAMB./2+(Gv./0.2e1)).*Hs.*LAMB).*F+Dd.*E.*Ld.*Cut.*Zcf1jx.*Zrf1jx.*rhocu.*(LAMB-LC+LC.*pcnt)).*A+E.*F.*Dd.*Hs.*Ld.*pcnt.*rcnt.*rhob.*Zcf1jx.*Zrf1jx.*LAMB).*alpha2+((eAe0+eBe0.*Prss).*F+E).*Zcf1jx.*((Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*rhocu.*Cut+Hs.*pcnt.*rcnt.*rhob.*Zrf1jx.*(Gv-LAMB+Ld)).*E.*A.*LAMB.*F.*Dd).*((eAe0+eBe0.*Prss).^3).*omega.^3-((eAe0+eBe0.*Prss).^2).*((((-2.*(Ld.*((-rhocu.*eBe0.*Prss./2+E.*Hs.*pcnt.*rcnt.*rhob+Zrf1jx.*Zcf1jx-rhocu.*eAe0./2).*LAMB-LC.*(pcnt-1).*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0)./2).*(eAe0+eBe0.*Prss).*Dd-E.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*LAMB.^2./2+((-LC.*eBe0.*rhocu.*(pcnt-1).*Prss./2-LC.*rhocu.*(pcnt-1).*eAe0./2+Zcf1jx.*((rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx)).*E-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)./2).*LAMB-LC.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)./2).*rhocu.*(Cut.^2)+(2.*((Zrf1jx.*Zcf1jx+E.*Hs.*pcnt.*rcnt.*rhob).*LAMB+LC.*Zrf1jx.*Zcf1jx.*(pcnt-1)).*Ld.*rhocu.*alpha1.*(eAe0+eBe0.*Prss).*Dd+((-Hs.*Prss.*eBe0.*pcnt.*rcnt.*rhob.*rhocu-Hs.*eAe0.*pcnt.*rcnt.*rhob.*rhocu+Zrf1jx.*Zcf1jx.*(2.*alpha1.*rhocu+Hs.*pcnt.*rcnt.*rhob)).*E-(eAe0+eBe0.*Prss).*Zrf1jx.*Zcf1jx.*rhocu).*LAMB.^2+((Gv.*Hs.*Prss.*eBe0.*pcnt.*rcnt.*rhob.*rhocu+Gv.*Hs.*eAe0.*pcnt.*rcnt.*rhob.*rhocu-Zcf1jx.*((-2.*alpha1.*(rhob.*rcnt.*Hs+LC.*Zrf1jx).*rhocu+Gv.*Hs.*rcnt.*rhob.*Zrf1jx).*pcnt+2.*LC.*Zrf1jx.*rhocu.*alpha1)).*E+(eAe0+eBe0.*Prss).*Zrf1jx.*Zcf1jx.*rhocu.*(-LC.*pcnt+LC+Gv)).*LAMB+Gv.*LC.*Zrf1jx.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Cut+E.*Hs.*pcnt.*rcnt.*rhob.*Zcf1jx.*Zrf1jx.*alpha1.*LAMB.*(Gv-LAMB)).*F+Cut.*E.*((Ld.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*Dd+Zcf1jx.*rhocu).*Cut+2.*Zrf1jx.*Zcf1jx.*((Dd.*Ld.*alpha1)-LAMB./2+(Gv./0.2e1))).*rhocu.*(LAMB-LC+LC.*pcnt)).*A+((Ld.*((Hs.*pcnt.*rcnt.*rhob.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*E+(eAe0+eBe0.*Prss).*Zrf1jx.*Zcf1jx.*rhocu).*LAMB+LC.*Zrf1jx.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Dd+(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB.*rhocu).*Cut+E.*Dd.*Hs.*Ld.*pcnt.*rcnt.*rhob.*Zcf1jx.*Zrf1jx.*alpha1.*LAMB).*F+Dd.*E.*Ld.*Cut.*Zcf1jx.*Zrf1jx.*rhocu.*(LAMB-LC+LC.*pcnt)).*(alpha2.^2)+((eAe0+eBe0.*Prss).*F+E).*(((-2.*(-E.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*LAMB.^2./2+((-LC.*eBe0.*rhocu.*(pcnt-1).*Prss./2-LC.*rhocu.*(pcnt-1).*eAe0./2+Zcf1jx.*((rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx)).*E-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)./2).*LAMB-LC.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)./2).*rhocu.*Dd.*(Cut.^2)+((((-Hs.*Prss.*eBe0.*pcnt.*rcnt.*rhob.*rhocu-Hs.*eAe0.*pcnt.*rcnt.*rhob.*rhocu+Zrf1jx.*Zcf1jx.*(2.*alpha1.*rhocu+Hs.*pcnt.*rcnt.*rhob)).*E-(eAe0+eBe0.*Prss).*Zrf1jx.*Zcf1jx.*rhocu).*LAMB.^2+((Hs.*eBe0.*pcnt.*rcnt.*rhob.*rhocu.*(Ld+Gv).*Prss+Hs.*pcnt.*rcnt.*rhob.*rhocu.*(Ld+Gv).*eAe0-Zcf1jx.*((-2.*alpha1.*(rhob.*rcnt.*Hs+LC.*Zrf1jx).*rhocu+Hs.*rcnt.*rhob.*Zrf1jx.*(Ld+Gv)).*pcnt+2.*LC.*Zrf1jx.*rhocu.*alpha1)).*E+(eAe0+eBe0.*Prss).*Zrf1jx.*Zcf1jx.*rhocu.*(-LC.*pcnt+LC+Ld+Gv)).*LAMB+LC.*Zrf1jx.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss).*(Ld+Gv)).*Dd+(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB.*rhocu).*Cut+Dd.*E.*Hs.*pcnt.*rcnt.*rhob.*Zcf1jx.*Zrf1jx.*alpha1.*LAMB.*(Gv-LAMB+Ld)).*F+(Cut.*rhocu+Zrf1jx.*(Gv-LAMB+Ld)).*Cut.*Zcf1jx.*E.*rhocu.*(LAMB-LC+LC.*pcnt).*Dd).*A+Cut.*(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB.*F.*rhocu.*Dd).*alpha2+Cut.*(((eAe0+eBe0.*Prss).*F+E).^2).*(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*A.*LAMB.*F.*rhocu.*Dd).*omega.^2+(-1.*1j).*(eAe0+eBe0.*Prss).*((((2.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx+E.*LAMB).*rhocu.*(LAMB-LC+LC.*pcnt).*(Cut.^2)+(-2.*Ld.*rhocu.*(LAMB+LC.*(pcnt-1)).*alpha1.*(eAe0+eBe0.*Prss).*Dd+((-Hs.*pcnt.*rcnt.*rhob-2.*alpha1.*rhocu).*E+rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB.^2+(((-2.*alpha1.*rhocu.*LC+rhob.*rcnt.*Hs.*Gv).*pcnt+2.*alpha1.*rhocu.*LC).*E-rhocu.*eBe0.*(-LC.*pcnt+LC+Gv).*Prss-rhocu.*(-LC.*pcnt+LC+Gv).*eAe0+(-LC.*pcnt.*Zrf1jx-2.*alpha1.*rhocu+Zrf1jx.*(LC+Gv)).*Zcf1jx).*LAMB-(Gv.*Prss.*eBe0.*rhocu+Gv.*eAe0.*rhocu-Zcf1jx.*(-2.*alpha1.*rhocu+Gv.*Zrf1jx)).*LC.*(pcnt-1)).*Cut-((Zrf1jx.*Zcf1jx+E.*Hs.*pcnt.*rcnt.*rhob).*LAMB+LC.*Zrf1jx.*Zcf1jx.*(pcnt-1)).*(Gv-LAMB).*alpha1).*(eAe0+eBe0.*Prss).*F+2.*E.*(rhocu.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx).*(Cut.^2)+(-Ld.*alpha1.*rhocu.*(eAe0+eBe0.*Prss).*Dd+(rhocu.*eBe0.*Prss./2+rhocu.*eAe0./2-Zrf1jx.*Zcf1jx./2).*LAMB-Gv.*Prss.*eBe0.*rhocu./2-Gv.*eAe0.*rhocu./2+Zcf1jx.*(-2.*alpha1.*rhocu+Gv.*Zrf1jx)./2).*Cut-Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB)./2).*(LAMB-LC+LC.*pcnt)).*A+((((Zrf1jx.*Zcf1jx-rhocu.*eAe0+E.*Hs.*pcnt.*rcnt.*rhob-rhocu.*eBe0.*Prss).*LAMB-LC.*(pcnt-1).*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx)).*Ld.*(eAe0+eBe0.*Prss).*Dd-E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB.^2+((-LC.*eBe0.*rhocu.*(pcnt-1).*Prss-LC.*rhocu.*(pcnt-1).*eAe0+Zcf1jx.*((rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx)).*E-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*LAMB-LC.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Cut-(((Zrf1jx.*Zcf1jx+E.*Hs.*pcnt.*rcnt.*rhob).*LAMB+LC.*Zrf1jx.*Zcf1jx.*(pcnt-1)).*Ld.*(eAe0+eBe0.*Prss).*Dd+(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB).*alpha1).*F-E.*(LAMB-LC+LC.*pcnt).*((Ld.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*Dd+Zcf1jx.*rhocu).*Cut+Dd.*Ld.*Zcf1jx.*Zrf1jx.*alpha1)).*(alpha2.^2)+((eAe0+eBe0.*Prss).*F+E).*(((2.*Dd.*rhocu.*(eAe0+eBe0.*Prss).*(LAMB-LC+LC.*pcnt).*(Zcf1jx+E.*LAMB).*(Cut.^2)+((((-Hs.*pcnt.*rcnt.*rhob-2.*alpha1.*rhocu).*E+rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB.^2+(((-2.*alpha1.*rhocu.*LC+rhob.*rcnt.*Hs.*(Ld+Gv)).*pcnt+2.*alpha1.*rhocu.*LC).*E-rhocu.*eBe0.*(-LC.*pcnt+LC+Ld+Gv).*Prss-rhocu.*(-LC.*pcnt+LC+Ld+Gv).*eAe0+(-LC.*pcnt.*Zrf1jx-2.*alpha1.*rhocu+Zrf1jx.*(LC+Ld+Gv)).*Zcf1jx).*LAMB-LC.*(rhocu.*eBe0.*(Ld+Gv).*Prss+rhocu.*(Ld+Gv).*eAe0-Zcf1jx.*(-2.*alpha1.*rhocu+Zrf1jx.*(Ld+Gv))).*(pcnt-1)).*(eAe0+eBe0.*Prss).*Dd-E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB.^2+((-LC.*eBe0.*rhocu.*(pcnt-1).*Prss-LC.*rhocu.*(pcnt-1).*eAe0+Zcf1jx.*((rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx)).*E-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*LAMB-LC.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Cut-alpha1.*((Gv-LAMB+Ld).*((Zrf1jx.*Zcf1jx+E.*Hs.*pcnt.*rcnt.*rhob).*LAMB+LC.*Zrf1jx.*Zcf1jx.*(pcnt-1)).*(eAe0+eBe0.*Prss).*Dd+(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB)).*F-E.*(-2.*rhocu.*Zcf1jx.*Dd.*(Cut.^2)+(((-rhocu.*eBe0.*Prss-rhocu.*eAe0+Zrf1jx.*Zcf1jx).*LAMB+rhocu.*eBe0.*(Ld+Gv).*Prss+rhocu.*(Ld+Gv).*eAe0-Zcf1jx.*(-2.*alpha1.*rhocu+Zrf1jx.*(Ld+Gv))).*Dd+Zcf1jx.*rhocu).*Cut+Dd.*Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB+Ld)).*(LAMB-LC+LC.*pcnt)).*A+(((-E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB.^2+((-LC.*eBe0.*rhocu.*(pcnt-1).*Prss-LC.*rhocu.*(pcnt-1).*eAe0+Zcf1jx.*((rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx)).*E-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*LAMB-LC.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Cut-(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB.*alpha1).*F-E.*rhocu.*Zcf1jx.*Cut.*(LAMB-LC+LC.*pcnt)).*Dd).*alpha2+(((eAe0+eBe0.*Prss).*F+E).^2).*A.*(((-E.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*LAMB.^2+((-LC.*eBe0.*rhocu.*(pcnt-1).*Prss-LC.*rhocu.*(pcnt-1).*eAe0+Zcf1jx.*((rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx)).*E-Zcf1jx.*rhocu.*(eAe0+eBe0.*Prss)).*LAMB-LC.*Zcf1jx.*rhocu.*(pcnt-1).*(eAe0+eBe0.*Prss)).*Cut-(Zrf1jx.*LAMB+(rhob.*rcnt.*Hs+LC.*Zrf1jx).*pcnt-LC.*Zrf1jx).*Zcf1jx.*E.*LAMB.*alpha1).*F-E.*rhocu.*Zcf1jx.*Cut.*(LAMB-LC+LC.*pcnt)).*Dd).*omega-(Cut-alpha1).*(((Gv-LAMB).*(eAe0+eBe0.*Prss).*((eAe0+eBe0.*Prss).*F+E).*A+((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx+E.*LAMB).*(eAe0+eBe0.*Prss).*F+E.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx)).*(alpha2.^2)+((((eAe0+eBe0.*Prss).*(Gv-LAMB+Ld).*Dd+Zcf1jx+E.*LAMB).*(eAe0+eBe0.*Prss).*F+((eAe0+eBe0.*Prss).*(Gv-LAMB+Ld).*Dd+Zcf1jx).*E).*A+((eAe0+eBe0.*Prss).*(Zcf1jx+E.*LAMB).*F+E.*Zcf1jx).*Dd).*((eAe0+eBe0.*Prss).*F+E).*alpha2+((eAe0+eBe0.*Prss).*(Zcf1jx+E.*LAMB).*F+E.*Zcf1jx).*A.*Dd.*(((eAe0+eBe0.*Prss).*F+E).^2)).*(LAMB-LC+LC.*pcnt)).*Hp+E.*LAMB.*F.*((-1.*1j).*(Cut.^2).*Zrf1jx.*Ld.*Zcf1jx.*A.*(alpha2.^2).*rhocu.^2.*((eAe0+eBe0.*Prss).^3).*Dd.*omega.^3-((((Ld.*(-2.*Zrf1jx.*Zcf1jx+rhocu.*eBe0.*Prss+rhocu.*eAe0).*Dd+Zcf1jx.*rhocu).*Cut+2.*Zrf1jx.*Zcf1jx.*((Dd.*Ld.*alpha1)-LAMB./2+(Gv./0.2e1))).*A+Ld.*Dd.*Zcf1jx.*Zrf1jx).*alpha2+(Cut.*rhocu+Zrf1jx.*(Gv-LAMB+Ld)).*((eAe0+eBe0.*Prss).*F+E).*Zcf1jx.*A.*Dd).*Cut.*alpha2.*rhocu.*((eAe0+eBe0.*Prss).^2).*omega.^2+(1j).*(((-2.*rhocu.*((Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx).*(Cut.^2)+(2.*Ld.*alpha1.*rhocu.*(eAe0+eBe0.*Prss).*Dd+(-rhocu.*eBe0.*Prss-rhocu.*eAe0+Zrf1jx.*Zcf1jx).*LAMB+Gv.*Prss.*eBe0.*rhocu+Gv.*eAe0.*rhocu-Zcf1jx.*(-2.*alpha1.*rhocu+Gv.*Zrf1jx)).*Cut+Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB)).*A+(Ld.*(rhocu.*eBe0.*Prss+rhocu.*eAe0-Zrf1jx.*Zcf1jx).*Dd+Zcf1jx.*rhocu).*Cut+Dd.*Ld.*Zcf1jx.*Zrf1jx.*alpha1).*(alpha2.^2)+((eAe0+eBe0.*Prss).*F+E).*((-2.*rhocu.*Zcf1jx.*Dd.*(Cut.^2)+(((-rhocu.*eBe0.*Prss-rhocu.*eAe0+Zrf1jx.*Zcf1jx).*LAMB+rhocu.*eBe0.*(Ld+Gv).*Prss+rhocu.*(Ld+Gv).*eAe0-Zcf1jx.*(-2.*alpha1.*rhocu+Zrf1jx.*(Ld+Gv))).*Dd+Zcf1jx.*rhocu).*Cut+Dd.*Zcf1jx.*Zrf1jx.*alpha1.*(Gv-LAMB+Ld)).*A+Cut.*rhocu.*Zcf1jx.*Dd).*alpha2+Cut.*(((eAe0+eBe0.*Prss).*F+E).^2).*Zcf1jx.*A.*rhocu.*Dd).*(eAe0+eBe0.*Prss).*omega-(Cut-alpha1).*(((eAe0+eBe0.*Prss).*(Gv-LAMB).*A+(Ld.*(eAe0+eBe0.*Prss).*Dd)+Zcf1jx).*(alpha2.^2)+((eAe0+eBe0.*Prss).*F+E).*(((eAe0+eBe0.*Prss).*(Gv-LAMB+Ld).*Dd+Zcf1jx).*A+Dd.*Zcf1jx).*alpha2+(((eAe0+eBe0.*Prss).*F+E).^2).*Zcf1jx.*A.*Dd)));

Y1c = -1./(omega.*imag(ZTTL));

end