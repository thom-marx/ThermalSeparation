within ThermalSeparation.Geometry.PackedColumn;
record MultipakI "Multipak-I"
//Mellapak 250Y; steel 1.4571
  extends Geometry(H=2, d=0.1, eps=0.654, rho_solid=7980*ones(n), c_solid= 500, a=370, sigma_crit=75e-3, zeta=0.1, d_char=0.025, theta=60, S=0.0075);

end MultipakI;
