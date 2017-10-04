function [TSFC, Ia, eta_th, eta_p, eta_o, m_dot] = prestazioni(ue, f, Qf, v0, T_req)
	Ia = (1 + f) * ue - v0;
	TSFC = f ./ Ia;
	eta_th = ((1 + f) * ue^2 - v0^2) / (2 * f * Qf);
	ni = v0 ./ ue;
	eta_p = 2 * ni / (1 + ni); 
	eta_o = eta_th * eta_p;
	m_dot = T_req ./ Ia;
end