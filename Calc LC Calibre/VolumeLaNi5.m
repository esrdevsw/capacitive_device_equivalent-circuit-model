% 
% %dados volume LaNi5Hx
% % EXP - PHYSICAL REVIEW B 64 184105		
% % Modelo 1  Physica B 403 (2008) 2372–2382		
% % Exp[27]  Physica B 403 (2008) 2372–2382 (Ono(1985))		
% % Teórico - PHYSICAL REVIEW B 64 184105
% 
% expPRB = [0.00317	86.73673
%             0.40719	87.19836
%             3.13053	96.13027
%             6.01183	108.82518
%             6.69848	108.82518];
%         
% Mod1PhB = [2	93.88
%             3	97.03
%             4	100.37];
%         
% Exp27ONO = [3	97.5
%             4	100
%             5	104
%             6	107.4];
%         
% TeoricPRB =[0.02645	85.98133
%             0.53188	87.70196
%             7.01604	108.44049]; 
%         
% % função BiDoseResp Nonlinear Curve Fit Origin        
%         
% HM = (0:0.013239999999996:7);
% 
% A1=	84.9124732953507;
% 	
% A2=	109.221575548341;
% 	
% LOGx01=	4.67879401097656;
% 	
% LOGx02=	1.34994199438351;
% 	
% h1=	0.669632969822042;
% 	
% h2=	0.653087764041975;
% 	
% ppar=	0.504761392638776;
%        
% span = A2 - A1;
% Section1 = span.*ppar./(1+10.^((LOGx01-HM).*h1));
% Section2 = span.* (1-ppar)./(1+10.^((LOGx02-HM).*h2));
%  
% BiDoseResp=A1 + Section1 +Section2;
%         
% figure    
% plot(expPRB(:,1),expPRB(:,2),'ko',...
%     Mod1PhB(:,1),Mod1PhB(:,2),'go',...
%     Exp27ONO(:,1),Exp27ONO(:,2),'bo',...
%     TeoricPRB(:,1),TeoricPRB(:,2),'mo',...
%     HM,BiDoseResp,'r--')
% grid on
% xlabel('x in LaNi_5H_x')
% ylabel('Volume [Angstrom ^3]')
% 
% legend('EXP - PHYSICAL REVIEW B 64 184105',...
% 'Modelo 1  Physica B 403 (2008) 2372–2382',	...	
% 'Exp[27]  Physica B 403 (2008) 2372–2382 (Ono(1985))'	,...
% 'Teórico - PHYSICAL REVIEW B 64 184105',...
% 'BiDoseResp','Location','best')


function DensTT = VolumeLaNi5(HM)
        
A1=	84.9124732953507;
	
A2=	109.221575548341;
	
LOGx01=	4.67879401097656;
	
LOGx02=	1.34994199438351;
	
h1=	0.669632969822042;
	
h2=	0.653087764041975;
	
ppar=	0.504761392638776;
       
span = A2 - A1;
Section1 = span.*ppar./(1+10.^((LOGx01-HM).*h1));
Section2 = span.* (1-ppar)./(1+10.^((LOGx02-HM).*h2));
 
BiDoseResp=A1 + Section1 +Section2;

Volcel=BiDoseResp.*1e-24; % cm3

% 	La	Ni	
ALa = 138.90547;%g/mol
ANi = 58.69340; %g/mol
% AH  = 1.008;

nLa = 1; %atomos
nNi = 5; %atomos
% nH = HM;

Nav = 6.02e23; %Numero de Avogrado

MsLa = ALa.*nLa;
MsNi = ANi.*nNi;
% MsH = AH.*nH./Nav;

% DensLa = MsLa./Volcel
% DensNi = MsNi./Volcel
% DensH  = MsH./Volcel

% DensTT = (MsLa + MsNi + MsH) ./Volcel

DensTT = (MsLa + MsNi ) ./(Nav.*Volcel);

% DensTT= 7.87;

end









