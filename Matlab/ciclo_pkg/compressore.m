function [p3, T3, s3, Lc] = compressore(p2, T2, beta_c, eta_ac)
	global gamma R Cp mdot
	p3 = beta_c * p2;
	T3p = T2 * beta_c ^ ((gamma - 1) / gamma);
	T3 = T2 + (T3p - T2) / eta_ac;
	s3 = sF(p3, T3, R, gamma);
	Lc = mdot * Cp * (T3 - T2);
end