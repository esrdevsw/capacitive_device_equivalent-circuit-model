function [xc0, pest, lambda, jacobian, residuals]=...
    meshLcRoTeta(Lcpar, Rhopar, P,C0,R0,m0,lmbd, frequ, HM)

global kp;

% global massa;
% global Prss

pico =1e-12;
kilo=1e3;

A = Lcpar;
B = Rhopar;

[AM, BM] = meshgrid(A,B) ;% generate 2D independent variables
n = size(AM);

% deltas 
q1 = (R0./kilo - R1PARarray( AM, BM, m0, P,lmbd, frequ, HM)./kilo).^2;
q2 = (C0./pico - C1PARarray( AM, BM, m0, P,lmbd, frequ, HM)./pico).^2;

minA=((min(min(min(q1)))));
minB=((min(min(min(q2)))));

min2A=(q1-minA);
min2B=(q2-minB);

maxA=((max(max(max(q1)))));
maxB=((max(max(max(q2)))));

Anorm=min2A./maxA;
Bnorm=min2B./maxB;

RCnorm = Anorm + Bnorm ;

Xminrc= min(min(RCnorm));
[XRC, YRC]= find(Xminrc==RCnorm);

LC1=mean(Lcpar(YRC));
rho0=mean(Rhopar(XRC));
disp('Zero malha')
disp([Xminrc LC1 rho0 ])

% pause(0.01)

figure(50);
clf
surfc(AM, BM, RCnorm,...'FaceLighting','phong',...
         'FaceColor','interp',...
         'EdgeColor','none');
    set(gca,'yscale','log')
    set(gca,'zscale','log')
    set(gca,'ZDir','reverse')
    colormap(flipud(colormap))
    % Create xlabel
    xlabel('Lc(mm)');
    % Create ylabel
    ylabel('\rho_b (\Omegam)');
    % Create zlabel
    zlabel('RC norm');
    shadowplot x
    shadowplot y
    
    view([-39 52]);
 
    
if rho0 > 10;
    RO1=rho0;
    kc=-20;
    while RO1 > 10
        kc=kc+1;
        RO1=rho0/10^kc ;
    end    
else
    RO1=rho0;
    kc=20;
    while (RO1 < 1) 
        kc=kc-1;
        RO1=rho0/10^kc ;
     end
end
    
kp = kc;
    
xc0=[LC1 RO1];

objfcn = @(xi)[R1PAR(xi, m0, P,lmbd, frequ, HM)/kilo - R0/kilo;...
               C1PAR(xi, m0, P,lmbd, frequ, HM)/pico - C0/pico]; 
    
opts = optimoptions(@lsqnonlin,...
        'Algorithm','trust-region-reflective',...
        'Display','off','MaxIter',1500 );%,...
%         'MaxFunEvals',500,...
%         'TolFun',1e-25, ...
%         'TolX',1e-25,...
%         'PlotFcns', {...
%         @optimplotx,... plots the current point
%         @optimplotfunccount... plots the function count
%         });
    
lb = xc0.*0.25; % Lower bound
ub = xc0.*1.95; % Upper bound

[pest,resnorm,residuals,exitflag,outputC, lambda, jacobian] = lsqnonlin(objfcn,xc0,lb,ub,opts);
 

end




