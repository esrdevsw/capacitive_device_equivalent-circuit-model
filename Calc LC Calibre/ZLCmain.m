% Multiplos arquivos Script de calculo das constantes da camara a partir de ensaios
% experimentias de camara vazia.

clear all
% close all
clc
format longg

menu =({'escolha entre os menus';...
    '<----------------------------------------->';
    '   menu     |     Function (data type)     ';...
    ' ----------------------------------------- ';...
    '     1 ==> carregar multiplos arquivos     ';...
    '     2 ==> carregar arquivo unico          ';...
    '     3 ==> calc. ponto CR(P,C,R,f e Lambda)';...
    '                                           ';...
    '                                           ';...    
    '                                           ';...
    ' ----------------------------------------- ';...
    '                                           ';...
    '                                           ';...
    '     9 ==> sair                            ';...
    '<----------------------------------------->'});

Data_Type = 0;

while Data_Type ~= 10
    
    disp (menu)
    prompt = 'Select = ';
    Data_Type = input(prompt);
    switch Data_Type
        case 1 
            %% calculo para multiplos arquivos
            delimiterIn = '\t';
            [filename,pathname,~] =  uigetfile({'CLC*.dat', 'All Files (*.*)'},'Pick a file','MultiSelect','on');
            if isequal(pathname,0)
                % No file Selected    
                disp('No file Selected!');
                return;
            elseif iscell(filename)
                % Multiple Files Selected        
                FilesCount=length(filename);                
                for loop1=1:FilesCount
                    sFileName=filename{loop1};        
                    if ~isequal(sFileName,0)            
                        sFullFileName=fullfile(pathname,sFileName);             
                        CTT = importdata(sFullFileName, delimiterIn);
                        nomeA=strsplit(sFileName,'_');
                        nomes=strsplit(nomeA{1},'lmbd');
                        nomeB=strsplit(nomeA{1},'lmbd');
                        Lambdas=str2double(nomeB{2});
                        %             [SET      HM       P        C         R         MZ        PHI     ]
                        completoMEAN=[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,11) (CTT(:,12)) CTT(:,13) CTT(:,14)];
           
                        frequ = mean(freq(completoMEAN(:,4), completoMEAN(:,6), completoMEAN(:,7)));
                        SET = completoMEAN(:,1);            
                        StartCiclos = find(SET==1); FimCiclos = length(SET);            
                        nciclos= length(StartCiclos);            
                        pontos=[StartCiclos; FimCiclos+1];           
                        A = SET;            
                        normA = max(A) - min(A);               % this is a vector            
                        normA = repmat(normA, [length(A) 1]);  % this makes it a matrix
                        Nciclos = A./normA ; % your normalized matrix

                        % CONFIRMAR NUMERO DE CICLOS DOS ENSAIOS
                        for ciclos = 1:nciclos
                            Svar = (pontos(ciclos):pontos(ciclos+1)-1);
                            Nciclos(Svar) = Nciclos(Svar)+ciclos;
                        end
            
                        figure('name',sFileName);

                        subplot(2,2,1)
                        plot(Nciclos,completoMEAN(:,4),'-*')
                        title('Capacidade')
                        ylabel('C [pF]');
                        xlabel('nº ciclos')
                        grid on

                        subplot(2,2,2)
                        plot(Nciclos,completoMEAN(:,5),'-*')
                        title('Resistividade')
                        ylabel('R [\Omega]');
                        xlabel('nº ciclos')
                        grid on

                        subplot(2,2,3)
                        plot(Nciclos,completoMEAN(:,6),'-*')
                        title('Módulo da impedância')
                        ylabel('|Z| [M\Omega]');
                        xlabel('nº ciclos')
                        grid on

                        subplot(2,2,4)
                        plot(Nciclos,completoMEAN(:,7),'-*')
                        title('Fase')
                        ylabel('\phi [º]');
                        xlabel('nº ciclos')
                        grid on
                    end
                    disp(sFullFileName)
                    disp({'Lambdas ',Lambdas, 'Nº Ciclos ',nciclos})
                    pause(1)
                end
    
                disp('<                                               >')
                disp('Entre com o intervalo dos ciclos para os arquivos')
                disp('<                                               >')
                %% intervalo dos ciclos a serem calculados
                for loop1=1:FilesCount
                    disp('<------------------------------------------------>')
                    sFileName=filename{loop1};
                    disp(sFileName)
                    prompt = 'Ciclo inicial = ';
                    inCLC(loop1) = input(prompt);        
                    prompt = 'Ciclo final = ';
                    FinCLC(loop1) = input(prompt); 
                    disp('<------------------------------------------------>')
                            menu2 =({'escolha entre os menus';...
                        '<----------------------------------------->';
                        '   menu     |     Function (data type)     ';...
                        ' ----------------------------------------- ';...
                        '     1 ==> Calc aquivo "medio"         ';...
                        '     2 ==> carregar aquivo "max"           ';...
                        '     3 ==> carregar aquivo "min"           ';...
                        '<----------------------------------------->'});

                    disp(menu2)
                    prompt = 'input = ';
                    TipoDados(loop1) = input(prompt);
                end
                %% calculo dos ciclos definidos
        
                for loop1=1:FilesCount
                    sFileName=filename{loop1};
                    if ~isequal(sFileName,0)
                        sFullFileName=fullfile(pathname,sFileName);            
                        CTT = importdata(sFullFileName, delimiterIn);  

                        nomeA=strsplit(sFileName,'_');
                        nomes=strsplit(nomeA{1},'lmbd');

                        nomeB=strsplit(nomeA{1},'lmbd');

                        lmbd=str2double(nomeB{2});

                        %              [SET      HM       P        C         R         MZ        PHI     ]
                        completoMEAN=[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,11) (CTT(:,12)) CTT(:,13) CTT(:,14)];
                        completoMIN =[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,16) (CTT(:,17)) CTT(:,18) CTT(:,19)]; 
                        completoMAX =[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,21) (CTT(:,22)) CTT(:,23) CTT(:,24)]; 

                        frequ = mean(freq(completoMEAN(:,4), completoMEAN(:,6), completoMEAN(:,7)));

                        SET = completoMEAN(:,1);

                        StartCiclos = find(SET==1); FimCiclos = length(SET);
                        nciclos= length(StartCiclos);
                        pontos=[StartCiclos; FimCiclos+1];

                        A = SET;
                        normA = max(A) - min(A);               % this is a vector
                        normA = repmat(normA, [length(A) 1]);  % this makes it a matrix

                        Nciclos = A./normA ; % your normalized matrix

                        % CONFIRMAR NUMERO DE CICLOS DOS ENSAIOS
                        for ciclos = 1:nciclos
                            Svar = (pontos(ciclos):pontos(ciclos+1)-1);
                            Nciclos(Svar) = Nciclos(Svar)+ciclos;
                        end
                        arquivo=TipoDados(loop1);
                        if arquivo == 1; datfile = completoMEAN; end
                        if arquivo == 2; datfile = completoMIN; end
                        if arquivo == 3; datfile = completoMAX; end

                        inicio = pontos(inCLC(loop1));
                        termino= pontos(FinCLC(loop1)+1)-1;

                        Theta=1;
                        [LC, RHO]=FunCalcCiclos(Nciclos,inicio,termino,datfile,sFileName)  
                    end
                end
                %% calculo para um unico arquivo
            else  %     para um unico arquivo
                sFileName=filename;
                sFullFileName=fullfile(pathname,sFileName);
            %     [~,sModelName,~] = fileparts(sModelName);
                CTT = importdata(sFullFileName, delimiterIn);

                nomeFile=strsplit(sFileName,'_');
                nomesLMB=strsplit(nomeFile{1},'lmbd');
                
%                 LAMBDA=str2double(nomesLMB{2});
                
%                 nomesMS=strsplit(nomeFile{2},'g');
                
%                 MASSA=str2double(nomesMS{1});            

                %              [SET      HM       P        C         R         MZ        PHI     ]
                completoMEAN=[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,11) (CTT(:,12)) CTT(:,13) CTT(:,14)];

%                 frequ = mean(freq(completoMEAN(:,4), completoMEAN(:,6), completoMEAN(:,7)));

                SET = completoMEAN(:,1);

                StartCiclos = find(SET==1); FimCiclos = length(SET);
                nciclos= length(StartCiclos);
                pontos=[StartCiclos; FimCiclos+1];

                A = SET;
                normA = max(A) - min(A);               % this is a vector
                normA = repmat(normA, [length(A) 1]);  % this makes it a matrix
                                                                       % of the same size as A
                Nciclos = A./normA ; % your normalized matrix

                % CONFIRMAR NUMERO DE CICLOS DOS ENSAIOS
                for ciclos = 1:nciclos
                    Svar = (pontos(ciclos):pontos(ciclos+1)-1);
                    Nciclos(Svar) = Nciclos(Svar)+ciclos;
                end

                figure('name',sFileName);

                subplot(2,2,1)
                plot(Nciclos,completoMEAN(:,4),'-*')
                title('Capacidade')
                ylabel('C [pF]');
                xlabel('nº ciclos')
                grid on

                subplot(2,2,2)
                plot(Nciclos,completoMEAN(:,5),'-*')
                title('Resistividade')
                ylabel('R [\Omega]');
                xlabel('nº ciclos')
                grid on

                subplot(2,2,3)
                plot(Nciclos,completoMEAN(:,6),'-*')
                title('Módulo da impedância')
                ylabel('|Z| [M\Omega]');
                xlabel('nº ciclos')
                grid on

                subplot(2,2,4)
                plot(Nciclos,completoMEAN(:,7),'-*')
                title('Fase')
                ylabel('\phi [º]');
                xlabel('nº ciclos')
                grid on

        
                disp(sFullFileName)
%                 disp({'Lambdas ',LAMBDA, 'Nº Ciclos ',nciclos})
                disp({ 'Nº Ciclos ',nciclos})
                pause(1)
    
                disp('<                                               >')
                disp('Entre com o intervalo dos ciclos para os arquivos')
                disp('<                                               >')
                %% intervalo dos ciclos a serem calculados
    
                disp('<------------------------------------------------>')
                sFileName=filename;
                disp(sFileName)
                prompt = 'Ciclo inicial = ';
                inCLC = input(prompt);        
                prompt = 'Ciclo final = ';
                FinCLC = input(prompt); 
                disp('<------------------------------------------------>')
                        menu2 =({'escolha entre os menus';...
                    '<----------------------------------------->';
                    '   menu     |     Function (data type)     ';...
                    ' ----------------------------------------- ';...
                    '     1 ==> Calc aquivo "medio"         ';...
                    '     2 ==> carregar aquivo "max"           ';...
                    '     3 ==> carregar aquivo "min"           ';...
                    '<----------------------------------------->'});

                disp(menu2)

                prompt = 'input = ';
                TipoDados = input(prompt);
                %% calculo dos ciclos definidos

                %              [SET      HM       P        C         R         MZ        PHI        Ciclo   massa    lambda]
                completoMEAN=[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,11) (CTT(:,12)) CTT(:,13) CTT(:,14) CTT(:,3) CTT(:,4) CTT(:,9)];
                completoMIN =[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,16) (CTT(:,17)) CTT(:,18) CTT(:,19) CTT(:,3)]; 
                completoMAX =[CTT(:,2) CTT(:,5) CTT(:,6) CTT(:,21) (CTT(:,22)) CTT(:,23) CTT(:,24) CTT(:,3)]; 
             
                arquivo=TipoDados;
                
                if arquivo == 1; datfile = completoMEAN; end
                if arquivo == 2; datfile = completoMIN; end
                if arquivo == 3; datfile = completoMAX; end

                inicio = pontos(inCLC);
                termino= pontos(FinCLC+1)-1;

%                 [LC, RHO]=FunCalcCiclos(Nciclos,inicio,termino,datfile,sFileName, MASSA, LAMBDA) 
                [LC, RHO, Theta]=FunCalcCiclos(Nciclos,inicio,termino,datfile,sFileName) ;
            end

        %% calculo ponto com CR
            case 3
                disp('<------------------------------------------------>')
                disp('<                                                >')
                disp('<Entre com os valores de P, C, R e freq. e Lambda>')
                disp('<                                                >')
                disp('<------------------------------------------------>')
                prompt = 'P = ';
                P = input(prompt);
                prompt = 'C = ';
                C0 = input(prompt);
                prompt = 'R = ';
                R0 = input(prompt);
                prompt = 'freq. = ';
                frequ = input(prompt);
                prompt = 'Lambda = ';
                lmbd = input(prompt);
                prompt= 'Massa = ';
                MSS = input(prompt);
                prompt= 'H/M = ';
                HM = input(prompt);
                                    
                [LC, RHO]=calcPontoCR(P,C0,R0,frequ,MSS, lmbd, HM);
                %% terminar programa
        case 9
            %terminar programa
            disp('terminar programa')
            button = questdlg('Ready to quit?', ...
            'Exit Dialog','Yes','No','No');
            switch button
                case 'Yes',
                    disp('Exiting MATLAB');
                    %Save variables to matlab.mat
                    quit cancel;
                    Data_Type = 10;
                case 'No',
                    disp('iniciar menus')
                    Data_Type = 0;
            end
            
        case 2
            %% calculo para um unico arquivo  
            
    end
end
