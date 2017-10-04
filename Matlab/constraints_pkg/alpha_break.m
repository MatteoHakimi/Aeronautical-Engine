function [a] = alpha_break(h, M0)
	global R gamma delta p_rif T_rif
	
	TR = 1;
	
	[~, ~, T_volo, p_volo] = tropos(h);
	theta_0 = T_volo / T_rif * (1 + delta * M0 ^ 2);
	delta_0 = p_volo / p_rif * (1 + delta * M0 ^ 2) ^ (gamma / (gamma - 1));
	
	if theta_0 <= TR
		a = delta_0 * (1 - 0.49 * sqrt(M0));
	else
		a = delta_0 * (1 - 0.49 * sqrt(M0) - 2 * (theta_0 - TR) + M0);
	end
end