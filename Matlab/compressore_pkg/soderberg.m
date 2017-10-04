
function [xi_rot, xi_stat] = soderberg(de_alfa, de_beta, h, chord, s, beta1, alfa2, W1, C2)
	nu = 3.5306e-5;
	
	% Rotore
	xi_rot = 0.04 + 0.06 * (de_beta / 100) ^ 2;
	xi_rot = (1 - xi_rot) * (0.975 + 0.075 * chord / h) - 1;
	
	D = 2 * h * s * cos(beta1) / (h + s * cos(beta1));
	Re = D * W1 / nu;

	xi_rot = (10^5 / Re) ^ (1/4) * xi_rot;

	% Statore
	xi_stat = 0.04 + 0.06 * (de_alfa / 100) ^ 2;
	xi_stat = (1 - xi_stat) * (0.975 + 0.075 * chord / h) - 1;
	
	D = 2 * h * s * cos(alfa2) / (h + s * cos(alfa2));
	Re = D * C2 / nu;

	xi_stat = (10^5 / Re) ^ (1/4) * xi_stat;
end