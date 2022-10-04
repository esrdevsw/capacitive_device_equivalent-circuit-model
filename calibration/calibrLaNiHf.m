clear all
clc

global frequ

pico =1e-12;
kilo=1e3;

format longg


% [FileName,PathName] = uigetfile('*.dat','Select the data file');

FileName=('Dados LaNi5Hf completo.dat');

CTT = importdata(FileName,'\t');   % '\t'-->reconhece valores separados por tab

SET = CTT(:,1);	
Msa = CTT(:,2);	% massa
Lc = CTT(:,3);  % Amostra mm
Pro = CTT(:,5); % porosidade
LBD = CTT(:,6); % lambda
Cap = CTT(:,7);	
Res = CTT(:,8);	
	
Prs = 2.14721;	% equivalente pressao para atm
frequ = 1001.73;

% figure
% plot3(LBD,Pro,Cap)
% grid on
% 
% figure
% plot3(LBD,Pro,Res)
% grid on

[CTTLmb,index] = sortrows(CTT, 6);

vetLBM=0;
nl=2.5;
for nlmb = 1:7
    
    valLmb= nl + 0.5;
    
    Ii = find(CTTLmb(:,6) == valLmb);
    
    lambdas(nlmb,:)=[valLmb min(Ii) max(Ii) max(Ii)-min(Ii) ];
    
    vetLBM=CTTLmb(min(Ii):max(Ii),:);
    
    [vetLBM,index] = sortrows(vetLBM, 5);
    
    eval([ 'vetLBM' num2str(nlmb) '=vetLBM;']);
    
%     figure(100)
%     hold on
%     plot(vetLBM(:,5),vetLBM(:,7),'o')
%     hold off
%     
%     figure(200)
%     hold on
%     plot(vetLBM(:,5),vetLBM(:,8),'o')
%     hold off

    nl = valLmb;
    
end

% Calculo da constante Bcnt para cada ponto sendo esta dependente de Lc,
% Lambda, Theta e tendo um Rho associado

Npt = length(CTT);
% for ii = 1:1    
for ii = 1:(Npt)
    
    Pctt = CTT(ii,:);
    lmbd=Pctt(6);
    R0= Pctt(8);
    C0= Pctt(7);
    
    VP=[];

    % limites dos parametros
    NBcn=250;
    Bcnmin=-5;
    Bcnmax=10;
    Bcnpar = logspace(Bcnmin,Bcnmax,NBcn);

    NRho=700;
    Rhomin=-13;
    Rhomax=11;
    Rhopar = (logspace(Rhomin,Rhomax,NRho));

    % paramentros zero
    loop = 0;
    emd=100;
    
    for jj = 1:50
                
        if jj>3 
            if emd > 5e-17
                disp('LIMITE DOS ERRO') 
                break 
            end
        end

%     while emd > 1e-15
        loop = loop+1;
        if loop > 10
            disp('LIMITE DOS LOOPS PARA MATRIZ')        
            break
        end
        % limites dos parametros
        if loop > 1
            NRho=100;
            NBcn=100;
        end
        if loop > 1
            disp ('entrou aqui')
            NBcn=NBcn+NBcn*loop*0.5;
            Bcnmin=Param(1)*0.5;
            Bcnmax=Param(1)*1.5;
            Bcnpar = linspace(Bcnmin,Bcnmax,NBcn);   

            Bcnpar = sort([Bcnpar Param(1)]);

            NRho=NRho+NRho*loop*0.5;
            Rhomin=(Param(2)*10^kc)*0.125;
            Rhomax=(Param(2)*10^kc)*2;
            Rhopar = linspace(Rhomin,Rhomax,NRho);

            Rhopar = sort([Rhopar Param(2)]);
        end

        [pest, lambda, jacobian, kc]=meshBcnRho(Bcnpar, Rhopar, Prs, Pctt);

        YcalLM   = [R1PAR(pest, Prs, Pctt);     C1PAR(pest, Prs, Pctt) ];

        VetPT = [R0; C0];

        % erro no ajuste do ponto
        e2 =sqrt((1./length(YcalLM)).*sum(((YcalLM-VetPT).^2)./(YcalLM.^2)));
        
%         DesvPad=0.2; % metade do erro da craveira(paquimetro)
%         graus=length(LCobj30)+length(ABC)-1;
%         chi2LC = (sum(((LCcalc(ABC) - Lc).^2)./DesvPad.^2))./graus;

        emd =  e2 ;

        vetPt=[e2 pest];

        figure(5000)        
        subplot(1,2,1)    
        plot(R0./kilo,R0./kilo,'-ko',R0./kilo,YcalLM(1)./kilo,'*-','markers',12)
%         xlim([R0*0.95, R0*1.05])
    %     ylim([R0*0.5, R0*1.5])
        grid on
        title('R (k\Omega)')
        subplot(1,2,2)
        plot(C0./pico,C0./pico,'-ko',C0./pico,YcalLM(2)./pico,'*-','markers',12)
%         xlim([C0*0.95, C0*1.05])
    %     ylim([min([C0,Ccal1])*0.95, min([C0,Ccal1])*1.5])
        grid on
        title('C (pF)')
        legend('Ponto','Lsqnon')

        VP =[VP; e2 pest]    

        VP = sortrows(VP);

        Param = VP(1,2:3);

        Bcn=Param(1);
        RHO=Param(2)*10^kc;


        disp('    ===========================')

        disp({'Bcn   =',num2str(Bcn);...
              'Rho  =',sprintf('%.4g', RHO);...
              'Erro =',sprintf('%.4g',VP(1,1))})

    end
    disp({'|==================','====================|';...
          '|     Bcn   =       ',[num2str(Bcn),' mm   '];...
          '|     Rho  =       ',[sprintf('%.4g',RHO),' Ohm m'];...
          '|==================','====================|'})
      
    SolVet(ii,:)=[Bcn RHO]
end

VETOR = [LBD Pro SolVet(:,1)];

BcnVet= SolVet(:,1);
RhoVet= SolVet(:,2);

 [fitresult, gof] = createFitBCNT(LBD, Pro, BcnVet)

