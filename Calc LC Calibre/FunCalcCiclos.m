function [LC, RHO, Theta]=FunCalcCiclos(Ncls,InCiclos,FimCiclos,datfile,sFileName) 
% global frequ
% global lmbd
% global Theta
% global expoente
kk=0;
for ii = InCiclos:FimCiclos
    
    kk=kk+1;
    
    Npt=ii;
    Indat = datfile(Npt,:);
%   [SET      HM       P        C         R         MZ        PHI        Ciclo   massa    lambda]
    HM = Indat(2);
    P = Indat(3);
    C0 = Indat(4);
    R0 = Indat(5);
    Z0 = Indat(6);
    F0 = Indat(7);
    ponto = Indat(8);
    MASSA = Indat(9);
    LAMBDA = Indat(10);

    frequ = 1001.73;%freq(C0, Z0, F0);
    
    [LCc, RHOc,DesvMalhaLC,DesvMalhaRO]=calcPontoCR(P,C0,R0,frequ,MASSA, LAMBDA, HM);
    
    
    
    [Y1r, Y1c, Tta] = RCfun(LCc, RHOc, MASSA, P ,LAMBDA, frequ, HM);
    
    LC(kk,:)   = LCc ;
    RHO(kk,:)  = RHOc ;
    Capft(kk,:)= Y1c;
    Resmft(kk,:)= Y1r;
    Theta(kk,:)= Tta;    
   
    
    Nciclos(kk,:) = ponto;	  
    Cap(kk,:) =  C0;	
    Res(kk,:) = R0;

        pause(0.1)
        
    figCLC=figure(111);   
    
    subplot(6,2,[1,3,5])
    plot(Nciclos,Capft,Nciclos,Cap,'o')
%     title('Capacidade')
    ylabel('C [F]')
    xlabel('nº ciclos')

    subplot(6,2,[7,9,11])
    plot(Nciclos,Resmft,Nciclos,Res,'o')
%     title('Resistividade')
    ylabel('R [\Omega]')
    xlabel('nº ciclos')
 
    subplot(6,2,[2,4])
    plot(Nciclos,LC,'-*')
%     title('Lc')
    ylabel('Lc [mm]');
%    xlabel('nº ciclos')
%    ylim([min(min(LC)),2*max(median(LC))])

    subplot(6,2,[6,8])
    plot(Nciclos,RHO,'-*')
%     title('Resistividade equivalente')
    ylabel('\rho_b [\Omegam]')
%     xlabel('nº ciclos')
    
    subplot(6,2,[10,12])
    plot(Nciclos,Theta,'-*')
%     title('Resistividade equivalente')
    ylabel('\Theta (porosidade)')
    xlabel('nº ciclos')
%     ylim([min(min(RHOc)),10*max(median(RHOc))])
 
    calcMean(kk,:) = [ P, LCc, RHOc ,C0,R0,Z0,F0, Tta];
    
    pause(0.1)
 
end

nomefigRc=['Ciclos_',sFileName,'.fig'];
saveas(figCLC,nomefigRc)
                            



        % open a file for writing
        fid = fopen((['Calc_param_',sFileName]), 'w');
        % print values in column order
        save((['Calc_param_',sFileName]),'calcMean','-ascii','-tabs','-append');
