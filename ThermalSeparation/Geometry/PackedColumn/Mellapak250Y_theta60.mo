within ThermalSeparation.Geometry.PackedColumn;
record Mellapak250Y_theta60 "Mellapak 250Y mit theta=60°"
//Mellapak 250Y; steel 1.4571
  extends Geometry(H=30, d=0.2, eps=0.98, rho_solid=7980*ones(n), c_solid= 500, a=250, sigma_crit=75e-3, zeta=0.1, d_char=0.025);//, theta=45);

end Mellapak250Y_theta60;
