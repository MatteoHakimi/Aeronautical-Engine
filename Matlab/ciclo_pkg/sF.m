function s = sF(p, T, R, gamma)
	Cp = R * gamma / (gamma - 1);
	s = Cp * log(T / 273.15) - R * log(p / 101325);
end