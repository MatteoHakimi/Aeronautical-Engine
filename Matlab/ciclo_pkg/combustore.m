function [p4, s4, f] = combustore(p3, T3, T4, eta_pb, eta_b, Qf)
	global gamma R Cp
	f = Cp * (T4 - T3) / (eta_b * Qf);
	p4 = eta_pb * p3;
	s4 = sF(p4, T4, R, gamma);
end