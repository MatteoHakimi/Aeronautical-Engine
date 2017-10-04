function [p2, T2, s2] = dinamica(p1, T1, eta_d, M)
	global R gamma
	delta = (gamma - 1) / 2;
	T2 = T1 * (1 + delta * M^2);
	T2p = T1 + eta_d * (T2 - T1);
	p2 = p1 * (T2p / T1) ^ (gamma / (gamma - 1));
	s2 = sF(p2, T2, R, gamma);
end