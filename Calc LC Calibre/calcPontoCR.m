function  [LC, RHO, DesvMalhaLC,DesvMalhaRO]=calcPontoCR(P0,C0,R0,frequ,MASSA, LAMBDA, HM)

global kp

VP=[];

% limites dos parametros
NLc=500;
Lcmin=0.01;
Lcmax=LAMBDA;
Lcpar = linspace(Lcmin,Lcmax,NLc);

NRho=700;
Rhomin=-2;
Rhomax=6;
Rhopar = (logspace(Rhomin,Rhomax,NRho));

% paramentros zero
loop = 0;
emd=100;

while emd > 1e-9
    loop = loop+1;
    if loop > 20
%         disp('LIMITE DOS LOOPS PARA MATRIZ')        
        break
    end
    % limites dos parametros
    if loop > 1
        NRho=100;
        NLc=200;
    end
    if loop > 1
%         disp ('entrou aqui')
        NLc=NLc+NLc*loop;
%         Lcmin=pest(1)*0.5;
%         Lcmax=pest(1)*1.5;
%         NLc=NLc+NLc*loop*0.025;
        Lcpar = linspace(Lcmin,Lcmax,NLc);   
            
        Lcpar = sort([Lcpar pest(1)]);
          
        NRho=NRho+NRho*loop*0.025;
        Rhomin=(pest(2)*10^kp)*0.125;
        Rhomax=(pest(2)*10^kp)*2;
        Rhopar = linspace(Rhomin,Rhomax,NRho);
            
        Rhopar = sort([Rhopar pest(2)]);
    end
    
    DesvMalhaLC=(Lcmax-Lcmin)/NLc;
    DesvMalhaRO=(Rhomax-Rhomin)/NRho;
    
    [xc0, pest, lambda, jacobian, residuals]=...
        meshLcRoTeta(Lcpar, Rhopar, P0,C0,R0,MASSA,LAMBDA, frequ, HM);
            
      
    Yxc0     = [R1PAR(xc0, MASSA, P0,LAMBDA, frequ, HM);      C1PAR(xc0, MASSA, P0,LAMBDA, frequ, HM) ];
 
    YcalLM   = [R1PAR(pest, MASSA, P0,LAMBDA, frequ, HM);     C1PAR(pest, MASSA, P0,LAMBDA, frequ, HM)    ];

    VetPT = [R0; C0];

    % erro no ajuste do ponto
    e0 =sqrt((1./length(Yxc0)).*sum(((Yxc0-VetPT).^2)./(Yxc0.^2)));
 
    e2 =sqrt((1./length(YcalLM)).*sum(((YcalLM-VetPT).^2)./(YcalLM.^2)));

    emd = mean([e0 e2 ]);

    vetPt=[e0 xc0; e2 pest];
        
    vetPt = sortrows(vetPt);
    
    Rcal1 = [R1PAR(xc0, MASSA, P0,LAMBDA, frequ, HM), R1PAR(pest, MASSA, P0,LAMBDA, frequ, HM)];
    Ccal1 = [C1PAR(xc0, MASSA, P0,LAMBDA, frequ, HM), C1PAR(pest, MASSA, P0,LAMBDA, frequ, HM)];
    
    figure(5000) 
    clf
    
    subplot(1,2,1)    
    plot(R0,R0,'-ko',R0,Rcal1,'*-','markers',12)
%     xlim([R0*0.95, R0*1.05])
%     ylim([R0*0.5, R0*1.5])
    grid on
    title('R')
    subplot(1,2,2)
    plot(C0,C0,'-ko',C0,Ccal1,'*-','markers',12)
%     xlim([C0*0.95, C0*1.05])
%     ylim([min([C0,Ccal1])*0.95, min([C0,Ccal1])*1.5])
    grid on
    title('C')
    legend('Ponto','mesh','Lsqnon')
    
%     pause(0.01)
    
    VP =[VP;e0 xc0; e2 pest]  ;  
    
    VP = sortrows(VP);
        
    Param = VP(1,2:3);
    
    LC=Param(1);
    RHO=Param(2)*10^kp;
        

    disp('    ===========================')
    
    disp({'Lc   =',num2str(LC);...
          'Rho  =',sprintf('%.4g', RHO);...
          'Erro =',sprintf('%.4g',VP(1,1))})

end
disp({'|==================','====================|';...
      '|     Lc   =       ',[num2str(LC),' mm   '];...
      '|     Rho  =       ',[sprintf('%.4g',RHO),' Ohm m'];...
      '|==================','====================|'})
end
    
