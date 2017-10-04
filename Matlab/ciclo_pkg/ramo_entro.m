function [s, T] = ramo_entro(p2, p1, T2, T1, s2, s1)
% 	global R gamma
% 	T = linspace(T1, T2, 100);
% 	p = linspace(p1, p2, 100);
% 	s = sF(p, T, R, gamma);
	global Cp R
	
	Tref = 293.15;
	pref = 101325;
% 	rhoref = pref / (R * Tref);
	
	rho2 = p2 / (R * T2);
	rho1 = p1 / (R * T1);
	
	n = log(p2 / p1) / log(rho2 / rho1);
	
	s = linspace(s2, s1, 100);
	
	T = (273.15 .* exp(s/Cp) .* (p1 / pref) .^ (R / Cp) .* (1 / T1) .^ (R * n / (Cp * (n - 1)))) .^ ((Cp * (n-1)) / (Cp * (n - 1) - R * n));
	
	s = flip(s);
	T = flip(T);
end