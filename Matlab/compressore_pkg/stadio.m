function [p2, p02, p3, p03, T2, T02, T3, T03, rho2, rho3, al_2_rm, be_1_rm, del_be, lam_rm, h, W1_rm, C2, M2, M3] = stadio(p1, p01, T1, T01, r_tip_1, del_h0_st, C_axial, omega, exp_n, m_dot)
	global R gamma delta C_P

	T03		= T01 + del_h0_st / C_P;
	T3		= T03 - C_axial ^ 2 / (2 * C_P);
	a3		= sqrt(gamma * R * T3);
	M3		= C_axial / a3;
	
	beta_st	= (T03 / T01) ^ (1 / exp_n);
	p03		= p01 * beta_st;
	p02		= p03;
	p3		= p03 / (1 + delta * M3 ^ 2) ^ (1 / exp_n);
	rho1		= p1 / (R * T1);
	A1		= m_dot / (rho1 * C_axial);
	
	chi		= sqrt(1 - A1 / (pi * r_tip_1 ^ 2));
	r_hub_1	= chi * r_tip_1;
	h		= r_tip_1 - r_hub_1;
	r_rm_1	= (r_tip_1 + r_hub_1) / 2;
	U_rm		= omega * r_rm_1;
	
	be_1_rm	= atan(U_rm / C_axial);
	be_2_rm	= atan(tan(be_1_rm) - del_h0_st / (U_rm * C_axial));
	lam_rm	= C_axial / (2 * U_rm) * (tan(be_1_rm) - tan(be_2_rm));
	
	W1_rm		= sqrt(C_axial ^ 2 + (C_axial * tan(be_1_rm)) ^ 2);
	T2		= T1 + del_h0_st/ C_P * lam_rm;
	a2		= sqrt(gamma * R * T2);
	al_2_rm	= atan(U_rm / C_axial - tan(be_2_rm));
	
	C2		= C_axial / cos(al_2_rm);
	M2		= C2 / a2;
	p2		= p02 / (1 + delta * M2 ^ 2) ^ (1/exp_n);
	T02		= T2 * (1 + delta * M2 ^ 2);
	
	rho2		= p2 / (R * T2);
	rho3		= p3 / (R * T3);
	
	del_be	= be_1_rm - be_2_rm;
end