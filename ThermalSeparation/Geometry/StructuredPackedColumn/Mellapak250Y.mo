within ThermalSeparation.Geometry.StructuredPackedColumn;
record Mellapak250Y "Mellapak 250Y for CCE"
//Mellapak 250Y; steel 1.4571
  extends Geometry(H=10, d=0.4, eps=0.96, rho_solid=7980*ones(n), c_solid= 500, a=250, sigma_crit=75e-3, zeta=0.1, d_char=0.025, theta=45);

end Mellapak250Y;
