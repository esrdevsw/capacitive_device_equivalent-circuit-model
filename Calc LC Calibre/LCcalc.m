function [LCobj30,ROobj] =  LCcalc(x)

global CTT
global P0
global m0
global C0
global R0
global lmbd

Lcexp=CTT(:,4);
Npt = length(Lcexp);

LCobj=zeros(Npt,1);
ii=0;
while ii<Npt   
    ii=ii+1;
    Indat = CTT(ii,5:11);
    
    P0 = Indat(1); m0 = Indat(2); C0 = Indat(3);  R0 = Indat(4); lmbd = Indat(7);   
    
    [LCc, RHOc,DesvMalhaLC,DesvMalhaRO]=calcPontoCR(x);   
    
    LCobj(ii,:)=LCc; 
    ROobj(ii,:)=RHOc; 
 
end

erroDesv = 0.2.*ones(size(Lcexp));
erroDesvMalhaLC = DesvMalhaLC.*ones(size(Lcexp));
erroDesvMalhaRO = DesvMalhaRO.*ones(size(Lcexp));

xpt=CTT(:,3);	
figure(1000)
% hold on
errorbar(xpt,Lcexp,erroDesv,'bo','markers',12)
grid on
title('LC fit')
xlabel('n ponto')
ylabel('LC [mm]')
legend('Ponto')
hold off

figure(1000)
hold on
errorbar(xpt,LCobj,erroDesvMalhaLC,'rs-','markers',6)
set(gca, 'YTick',(0:0.4:7)) 
grid on
title('LC fit')
xlabel('n ponto')
ylabel('LC [mm]')
legend('Ponto','Lc Calc','Location','NorthWest')
hold off


Nciclos = CTT(:,3);	
Lc = CTT(:,4);	
P = CTT(:,5);	
Msa = CTT(:,6);	
Cap = CTT(:,7);	
Res = CTT(:,8);	
LAMBDA = CTT(:,11);


[ResC, Ttar] = R1PARarray(x, LCobj, ROobj, Msa, P, LAMBDA);
    
[CapC, Ttac] = C1PARarray(x, LCobj, ROobj, Msa, P, LAMBDA);

       
    figCLC=figure(111);   
    subplot(6,2,[1,3,5])
    plot(Nciclos,CapC,Nciclos,Cap,'o')
%     title('Capacidade')
    ylabel('C [pF]')
    xlabel('nº ciclos')

    subplot(6,2,[7,9,11])
    plot(Nciclos,ResC,Nciclos,Res,'o')
%     title('Resistividade')
    ylabel('R [\Omega]')
    xlabel('nº ciclos')
 
    subplot(6,2,[2,4])
    plot(Nciclos,LCobj,'-*' ,Nciclos,Lc,'o')
%     title('Lc')
    ylabel('Lc [mm]');
%    xlabel('nº ciclos')
%    ylim([min(min(LC)),2*max(median(LC))])

    subplot(6,2,[6,8])
    errorbar(Nciclos,ROobj,erroDesvMalhaRO,'-*')
%     title('Resistividade equivalente')
    ylabel('\rho_b [\Omegam]')
%     xlabel('nº ciclos')
    
    subplot(6,2,[10,12])
    plot(Nciclos,Ttac,'-*')
%     title('Resistividade equivalente')
    ylabel('\Theta (porosidade)')
    xlabel('nº ciclos')
%     ylim([min(min(RHOc)),10*max(median(RHOc))])


LCobj30=LCobj;


