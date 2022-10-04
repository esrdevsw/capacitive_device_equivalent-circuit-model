% constantes globais
function [rcnt, pcnt, cbcnt]=  FUNCNTs(Tta,Bcnt)

% rcnt = Acnt.^(1-Tta);
% pcnt = Bcnt./Tta;
% cbcnt= Ccnt.*Tta.^2;
%% -----------------------------------------------------------
% rcnt = Acnt.*((1-Tta)) ;
% cbcnt= Ccnt.*(Tta);
% pcnt = Bcnt.^(cbcnt./rcnt);
%% -----------------------------------------------------------
% rcnt = 1;
% 
% cbcnt= Ccnt.*(Tta);
% 
% pcnt = Bcnt.^(cbcnt./((1-Tta)));
%% -----------------------------------------------------------
% rcnt = 1;
% 
% cbcnt= Ccnt.*(Tta);
% 
% pcnt = Bcnt.^(Tta./(Ccnt.*(1-Tta)));
%% -----------------------------------------------------------
% rcnt = 1;
% 
% cbcnt= Ccnt.*(Tta);
% 
% pcnt = Bcnt.^(((1-Tta)./(Tta)));
%% -----------------------------------------------------------
rcnt = 1;

cbcnt= 1;

pcnt = Tta.*Bcnt;

