within ThermalSeparation.Geometry.StructuredPackedColumn;
record Rombopak_6M "Rombopak-6M"
//Mellapak 250Y; steel 1.4571
  extends Geometry(H=1, d=0.1, eps=0.95, rho_solid=7980*ones(n), c_solid= 500, a=230, sigma_crit=75e-3, zeta=0.1, d_char=0.025, theta=60);

end Rombopak_6M;
