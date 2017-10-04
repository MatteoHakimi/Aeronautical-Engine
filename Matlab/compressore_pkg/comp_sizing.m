function [p1_vec, T1_vec, p2_vec, T2_vec, p3_vec, T3_vec, M2, M3, lam_vec, num_stadi, num_pale, xi_rot, xi_stat, h_vec, r_tip_in] = comp_sizing(p_in, T_in, M_in, M_tip_in, beta_tot, chi, AR, DF, m_dot, eta_p, lam_rm_in)
	global R gamma delta C_P
	
	r2d = @(x) x * 180 / pi;

	n_poli = gamma * eta_p / (1 + gamma * (eta_p - 1));
	exp_p = (gamma - 1) / gamma;
	exp_n = (n_poli - 1) / n_poli;
	
	a_in		= sqrt(gamma * R * T_in);
	rho_in	= p_in / (R * T_in);
	
	C_axial	= M_in * a_in;
	T0_in		= T_in * (1 + delta * M_in ^ 2);
	p0_in		= p_in * (1 + delta * M_in ^ 2) ^ (1/exp_p);
	T0_out	= T0_in * beta_tot ^ exp_n;
	del_h0_tot	= C_P * (T0_out - T0_in);
	
	A_in		= m_dot / (rho_in * C_axial);
	r_tip_in	= sqrt(A_in / (pi * (1 - chi ^ 2)));
	r_hub_in	= chi * r_tip_in;
	r_rm_in	= (r_tip_in +  r_hub_in) / 2;
	
	W1_tip_in	= M_tip_in * a_in;
	U_tip_in	= sqrt(W1_tip_in ^ 2 - C_axial ^ 2);
	omega		= U_tip_in / r_tip_in;
	U_rm_in	= omega * r_rm_in;
	U_hub_in	= omega * r_hub_in;
	
	W1_rm_in	= sqrt(U_rm_in ^ 2 + C_axial ^ 2);
	W1_hub_in	= sqrt(U_hub_in ^ 2 + C_axial ^ 2);
	
	be_1_rm_in	= acos(C_axial / W1_rm_in);
	be_1_tip_in	= acos(C_axial / W1_tip_in);
	be_1_hub_in	= acos(C_axial / W1_hub_in);
	
	check_lam	= 1 - (r_hub_in / r_rm_in) ^ 2;
	if (lam_rm_in < check_lam)
		lam_rm_in = check_lam;
		correz_lam = true;	% Altrimenti lambda negativo all'hub in
	else
		correz_lam = false;
	end
	
	be_2_rm_in	= atan(2 * U_rm_in * lam_rm_in / C_axial - tan(be_1_rm_in));
	del_h0_st	= U_rm_in * C_axial * (tan(be_1_rm_in) - tan(be_2_rm_in));
	
	num_stadi	= ceil(del_h0_tot / del_h0_st);
	del_h0_st	= del_h0_tot / num_stadi;
	
	be_2_rm_in	= atan(tan(be_1_rm_in) - del_h0_st / (U_rm_in * C_axial));
	
	lam_rm_in	= C_axial / (2 * U_rm_in) * (tan(be_1_rm_in) - tan(be_2_rm_in));
	
	if (lam_rm_in < check_lam)
		num_stadi = num_stadi + 1;
		correz_stadi = true;	% Altrimenti lambda negativo all'hub in
	end
	
	del_h0_st	= del_h0_tot / num_stadi;
	be_2_rm_in	= atan(tan(be_1_rm_in) - del_h0_st / (U_rm_in * C_axial));
	
	lam_rm_in	= C_axial / (2 * U_rm_in) * (tan(be_1_rm_in) - tan(be_2_rm_in));
	lam_hub_in	= 1 - (1 - lam_rm_in) / (r_hub_in / r_rm_in) ^ 2;
	lam_tip_in	= 1 - (1 - lam_rm_in) / (r_tip_in / r_rm_in) ^ 2;
	
	be_2_tip_in	=  atan(2 * U_tip_in * lam_tip_in / C_axial - tan(be_1_tip_in));
	be_2_hub_in	=  atan(2 * U_hub_in * lam_hub_in / C_axial - tan(be_1_hub_in));
	
	al_2_rm_in	= atan(U_rm_in / C_axial - tan(be_2_rm_in));
	al_2_hub_in	= atan(U_hub_in / C_axial - tan(be_2_hub_in));
	al_2_tip_in	= atan(U_tip_in / C_axial - tan(be_2_tip_in));
	
	W2_u_rm_in	= C_axial * tan(be_2_rm_in);
	W2_u_hub_in	= C_axial * tan(be_2_hub_in);
	W2_u_tip_in	= C_axial * tan(be_2_tip_in);
	
	C2_u_rm_in	= C_axial * tan(al_2_rm_in);
	C2_u_hub_in	= C_axial * tan(al_2_hub_in);
	C2_u_tip_in	= C_axial * tan(al_2_tip_in);
	
	W2_rm_in	= sqrt(W2_u_rm_in ^ 2 + C_axial ^ 2);
	W2_hub_in	= sqrt(W2_u_hub_in ^ 2 + C_axial ^ 2);
	W2_tip_in	= sqrt(W2_u_tip_in ^ 2 + C_axial ^ 2);
	
	C2_rm_in	= sqrt(C2_u_rm_in ^ 2 + C_axial ^ 2);
	C2_hub_in	= sqrt(C2_u_hub_in ^ 2 + C_axial ^ 2);
	C2_tip_in	= sqrt(C2_u_tip_in ^ 2 + C_axial ^ 2);
	
	si_rm_in_r	= C2_u_rm_in / (2 * (DF - 1) * W1_rm_in + W2_rm_in);
	si_hub_in_r	= C2_u_hub_in / (2 * (DF - 1) * W1_hub_in + W2_hub_in);
	si_tip_in_r	= C2_u_tip_in / (2 * (DF - 1) * W1_tip_in + W2_tip_in);
	
	si_rm_in_s	= C2_u_rm_in / (2 * (DF - 1) * C_axial + C2_rm_in);
	si_hub_in_s	= C2_u_hub_in / (2 * (DF - 1) * C_axial + C2_hub_in);
	si_tip_in_s	= C2_u_tip_in / (2 * (DF - 1) * C_axial + C2_tip_in);
	
	h_annul_in	= r_tip_in - r_hub_in;
	c_in		= AR * h_annul_in;
	
	si_vec	= [si_rm_in_r, si_rm_in_s, si_hub_in_r, si_hub_in_s, si_tip_in_r, si_tip_in_s];
	r_vec		= [r_rm_in, r_rm_in, r_hub_in, r_hub_in, r_tip_in, r_tip_in];
	
	[si_max, i_si_max] = max(si_vec);
	r_design	= r_vec(i_si_max);
	s_min		= c_in / si_max;
	
	num_pale	= ceil(2 * pi * r_design / s_min);
	
	% Vettori stadio
	T1_vec = zeros(1, num_stadi);		T1_vec(1) = T_in;
	T2_vec = zeros(1, num_stadi);	
	T3_vec = zeros(1, num_stadi);
	p1_vec = zeros(1, num_stadi);		p1_vec(1) = p_in;
	p2_vec = zeros(1, num_stadi);
	p3_vec = zeros(1, num_stadi);
	rho1_vec = zeros(1, num_stadi);	rho1_vec(1) = rho_in;
	rho2_vec = zeros(1, num_stadi);
	rho3_vec = zeros(1, num_stadi);
	T01_vec = zeros(1, num_stadi);	T01_vec(1) = T0_in;
	T02_vec = zeros(1, num_stadi);
	T03_vec = zeros(1, num_stadi);
	p01_vec = zeros(1, num_stadi);	p01_vec(1) = p0_in;
	p02_vec = zeros(1, num_stadi);
	p03_vec = zeros(1, num_stadi);
	be_1_vec = zeros(1, num_stadi);	be_1_vec(1) = be_2_rm_in;
	del_be_vec = zeros(1, num_stadi);	del_be_vec(1) = be_1_rm_in - be_2_rm_in;
	lam_vec = zeros(1, num_stadi);	lam_vec(1) = lam_rm_in;
	del_al_vec = zeros(1, num_stadi);	del_al_vec(1) = al_2_rm_in;
	h_vec = zeros(1, num_stadi);		h_vec(1) = h_annul_in;
	xi_stat = zeros(1, num_stadi);
	xi_rot = zeros(1, num_stadi);
	M2 = zeros(1, num_stadi);
	M3 = zeros(1, num_stadi);

	
	for i = 1 : num_stadi
		[p2_vec(i), p02_vec(i), p3_vec(i), p03_vec(i), T2_vec(i), T02_vec(i), T3_vec(i), T03_vec(i), rho2_vec(i), rho3_vec(i), del_al_vec(i), be_1_vec(i), del_be_vec(i), lam_vec(i), h_vec(i), W1, C2, M2(i), M3(i)] = ...
			stadio(p1_vec(i), p01_vec(i), T1_vec(i), T01_vec(i), r_tip_in, del_h0_st, C_axial, omega, exp_n, m_dot);
		T1_vec(i + 1) = T3_vec(i);
		p1_vec(i + 1) = p3_vec(i);
		rho1_vec(i + 1) = rho3_vec(i);
		T01_vec(i + 1) = T03_vec(i);
		p01_vec(i + 1) = p03_vec(i);


		[xi_rot(i), xi_stat(i)] = soderberg(del_al_vec(i), del_be_vec(i), h_vec(i), c_in, s_min, be_1_vec(i), del_al_vec(i), W1, C2);
	end
end
		