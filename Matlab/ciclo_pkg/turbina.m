function [p5, T5, s5] = turbina(p4, T4, Lc, eta_at, eta_mt, eta_mc)
	global R gamma Cp mdot
	Lt = Lc / (eta_mt * eta_mc);
	T5 = T4 - Lt / (Cp * mdot);
	T5p = T4 - (T4 - T5) / eta_at;
	p5 = p4 * (T5p / T4) ^ (gamma / (gamma - 1));
	s5 = sF(p5, T5, R, gamma);
end