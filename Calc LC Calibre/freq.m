function freqCalc = freq(Cap, ModZ, fase)
% calculo da impedancia complexa
modZ = ModZ; 
Phi = degtorad(fase);
Zc = modZ.*exp(1j.*Phi);

Iz = (Zc-conj(Zc))./(2*1j);

ff = (-1)./(Iz.*Cap.*2.*pi);

Npnt = length(ff);				  % number of data points

t = (1:Npnt)';				  % independent variable

freqCalc = polyfit(t,ff,0);