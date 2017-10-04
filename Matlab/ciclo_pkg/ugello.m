function [T9, ue, s9] = ugello(p5, T5, pa, eta_n)
	global R gamma Cp
	T9 = T5 * (1 - eta_n * (1 - (pa ./ p5) ^ ((gamma - 1) / gamma)));
	ue = sqrt(2 * Cp * (T5 - T9));
	s9 = sF(pa, T9, R, gamma);
end