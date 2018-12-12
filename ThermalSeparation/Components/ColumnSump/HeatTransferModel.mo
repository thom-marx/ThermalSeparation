within ThermalSeparation.Components.ColumnSump;
package HeatTransferModel
  model BaseHeatTransferModel
  end BaseHeatTransferModel;

  model dT_m "heat transfer calc. with dT_m; alpha constant"

  extends
      ThermalSeparation.Components.ColumnSump.HeatTransferModel.BaseHeatTransferModel;

  output SI.HeatFlowRate Q_kon;

  /* fluid properties */
  parameter Integer ns=2;
  input SI.Pressure p_sump;
  input SI.MoleFraction x_l[ns];
  input SI.MassFlowRate mdot_l_in;

  /* steam inflow */
  input SI.Pressure p_s;
  input SI.VolumeFlowRate vdot_v_in;
  SI.Density rho_v_in(start=0.6)=Water.dewDensity(sat);
  SI.Temperature T_s_in=T_s_out;
  SI.SpecificEnthalpy h_v_in=Water.dewEnthalpy(sat);
  SI.SpecificEnthalpy u_v_in=h_v_in-p_s/rho_v_in;

  /* steam outflow */
  SI.VolumeFlowRate vdot_v_out(start=0.5);
  SI.VolumeFlowRate vdot_l_out;
  SI.Density rho_v_out(start=0.6)=Water.dewDensity(sat);
  SI.Density rho_l_out(start=1000)=Water.bubbleDensity(sat);
  SI.Temperature T_s_out(start=T_s_in);
  SI.SpecificEnthalpy h_l_out(start=1e3)=Water.bubbleEnthalpy(sat);
  SI.SpecificEnthalpy h_v_out(start=1e4)=Water.dewEnthalpy(sat);
  SI.SpecificEnthalpy u_l(start=1e3)=h_l_out-p_s/rho_l_out;
  SI.SpecificEnthalpy u_v(start=1e4)=h_v_out-p_s/rho_v_out;

  /* heat transfer */
  SI.TemperatureDifference dT_m(start=1);
  input SI.Temperature T_f_in;
  input SI.Temperature T_f_out;
  SI.CoefficientOfHeatTransfer k;
  parameter SI.CoefficientOfHeatTransfer alpha_i=1e3;
  parameter SI.CoefficientOfHeatTransfer alpha_a=1e4;

  /* geometry */
  parameter SI.ThermalConductivity lambda_wall=15;
  parameter SI.Length L=6; //tube length
  parameter SI.Length d_a=35/1000; //outer tube diameter
  parameter SI.Length d_i=31/1000; //inner tube diameter
  parameter Real n=160; //number of tubes inside the heat exchanger
  parameter SI.Length s=(d_a-d_i)/2; //wall thickness of tube
  SI.Length d_m;
  SI.Area A_a; // inner surface area of the tubes
  SI.Area A_i; // outer surface area of the tubes
  SI.Area A_m; // mean value of surface area
  SI.Length h_w=L/20; //weir height
  SI.Length h_lw=0.7*((pi/4*d_i^2)*L*n); //weir length
    import Modelica.Constants.pi;

  SI.Volume HU_l;
  SI.Volume HU_v;

  SI.Length h_liq;

  /* medium */

  replaceable package Water = 
  Modelica.Media.Water.StandardWater 
  constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium annotation(choicesAllMatching=true);

  Water.SaturationProperties sat
      "State vector to compute saturation properties";

  equation
  /* geometry*/
   A_a= pi*d_a*L*n;
   A_i= pi*d_i*L*n;
   A_m=(A_a-A_i)/ln(A_a/A_i);

   d_m=(d_a-d_i)/ln(d_a/d_i);

  (pi/4*d_i^2)*L*n=HU_l+HU_v;
  h_liq=HU_l/((pi/4*d_i^2)*L*n);

  /* outflow */
  vdot_l_out=1.848*h_lw*(abs(h_liq-h_w))^1.5; //outflow equation

  /* mass balance */
  der(HU_l*rho_l_out+HU_v*rho_v_out)=vdot_v_in*rho_v_in-vdot_v_out*rho_v_out-vdot_l_out*rho_l_out;

  /* energy balance */
  der(u_l*rho_l_out*HU_l+u_v*rho_v_out*HU_v)=vdot_v_in*h_v_in*rho_v_in-vdot_v_out*h_v_out*rho_v_out-vdot_l_out*h_l_out*rho_l_out-Q_kon;

  /* saturation temperature */
  T_s_out=Water.saturationTemperature(p_s);
  sat.Tsat=T_s_out;
  sat.psat=p_s;
  //der(p_s)/der(T_s_out)=2257e3/T_s_out*(HU_v-HU_l); //clausius clapeyron

  /* heat transfer */
  Q_kon=k*A_m*dT_m;

  (T_s_in-T_f_out)/(T_s_out-T_f_in)=exp(((T_s_in-T_f_out)-(T_s_out-T_f_in))/dT_m); //mean temperature difference for countercurrent flow

  k=1/((1/alpha_i)*(d_a/d_i)+(s/lambda_wall)*(d_a/d_m)+1/alpha_a);

  initial equation
  //T_out=273+120;
  //HU_l=0.05*A*H;
  h_liq=h_w;
  //vdot_v_out=vdot_v_in;

  assert(Q_kon> 0, "Q_kon< 0; heat is flowing in the wrong direction");

    annotation (Diagram(graphics), Icon(graphics));
  end dT_m;

  model dT_m_2 "heat transfer calc. with dT_m; alpha not constant"

  extends
      ThermalSeparation.Components.ColumnSump.HeatTransferModel.BaseHeatTransferModel;

  output SI.HeatFlowRate Q_kon;

  /* fluid properties */
  parameter Integer ns=2;
  input SI.Pressure p_sump;
  input SI.MoleFraction x_l[ns];
  input SI.MassFlowRate mdot_l_in;

  /* steam inflow */
  input SI.Pressure p_s;
  input SI.VolumeFlowRate vdot_v_in;
  SI.Density rho_v_in(start=0.6)=Water.dewDensity(sat);
  SI.Temperature T_s_in=T_s_out;
  SI.SpecificEnthalpy h_v_in=Water.dewEnthalpy(sat);
  SI.SpecificEnthalpy u_v_in=h_v_in-p_s/rho_v_in;

  /* steam outflow */
  SI.VolumeFlowRate vdot_v_out(start=0.5);
  SI.VolumeFlowRate vdot_l_out;
  SI.Density rho_v_out(start=0.6)=Water.dewDensity(sat);
  SI.Density rho_l_out(start=1000)=Water.bubbleDensity(sat);
  SI.Temperature T_s_out(start=T_s_in);
  SI.SpecificEnthalpy h_l_out(start=1e3)=Water.bubbleEnthalpy(sat);
  SI.SpecificEnthalpy h_v_out(start=1e4)=Water.dewEnthalpy(sat);
  SI.SpecificEnthalpy u_l(start=1e3)=h_l_out-p_s/rho_l_out;
  SI.SpecificEnthalpy u_v(start=1e4)=h_v_out-p_s/rho_v_out;

  /* heat transfer */
  SI.TemperatureDifference dT_m(start=1);
  input SI.Temperature T_f_in;
  input SI.Temperature T_f_out;
  SI.CoefficientOfHeatTransfer k;
  SI.SpecificEnthalpy h_evap(start=1e6)=h_v_out-h_l_out;

  /* geometry */
  parameter SI.ThermalConductivity lambda_wall=15;
  parameter SI.Length L=6; //tube length
  parameter SI.Length d_a=35/1000; //outer tube diameter
  parameter SI.Length d_i=31/1000; //inner tube diameter
  parameter Real n=160; //number of tubes inside the heat exchanger
  parameter SI.Length s=(d_a-d_i)/2; //wall thickness of tube
  SI.Length d_m;
  SI.Area A_a; // inner surface area of the tubes
  SI.Area A_i; // outer surface area of the tubes
  SI.Area A_m; // mean value of surface area
  SI.Length h_w=L/20; //weir height
  SI.Length h_lw=0.7*((pi/4*d_i^2)*L*n); //weir length
    import Modelica.Constants.pi;

  SI.Volume HU_l;
  SI.Volume HU_v;

  SI.Length h_liq;

  constant Real g=Modelica.Constants.g_n;

  /* medium */

  replaceable package Water = 
  Modelica.Media.Water.StandardWater 
  constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium annotation(choicesAllMatching=true);

  Water.SaturationProperties sat
      "State vector to compute saturation properties";

  SI.SpecificHeatCapacity cp_l_out = 4187;
  SI.SpecificHeatCapacity cp_v_out = 2027;

  /* for calculation ov average coefficient of heat transfer */

  Real sigma; //coefficient of drag without condensation (1-phase)
  Real sigma_lam;
  Real sigma_turb;
  Real sigma_star; //coefficient of drag with condensation (2-phase)
  SI.ShearStress tau_s; //shear stress of steam
  SI.ShearStress tau_s_star; //normalized shear stress
  SI.Velocity v(start=0.1); //velocity of steam
  Real E; //         amplification factor
  Real psi;
  SI.CoefficientOfHeatTransfer alpha_i(start=1e2);// valid for laminar and turbulent flow; VDI-Wärmeatlas
  SI.CoefficientOfHeatTransfer alpha_a(start=1e4);//for bulk boiling in pressure vessels (natural convection); Stephan, Baehr, "Wärme- und Stoffübertragung"
  Real Nu_turb_total;//Nusselt number; valid for laminar and turbulent flow
  Real Nu_lam;
  Real Nu_lam_1;
  Real Nu_lam_2;
  Real Nu_turb;
  Real Re_f(start=1); //Reynolds number of the condensate film
  Real Re_s; //Reynolds number of steam
  SI.Length L_char; //characteristic length
  SI.KinematicViscosity ny_water;
  SI.KinematicViscosity ny_steam;
  SI.DynamicViscosity eta_water=4.65e-4;
  SI.DynamicViscosity eta_steam=1.295e-5;
  Real Tau;
  Real Ga; //Galilei number
  Real Ph; //phase change number
  Real Pr; //Prandtl number
  Real a;
  Real b;
  Real c;
  Real d;
  Real e;
  Real f;
  Real m; //tabulated value - VDI Wärmeatlas
  Real f_well;
  SI.MassFlowRate mdot_cond;
  SI.ThermalConductivity lambda_water=0.66;
  SI.Temperature T_wall(start=273+100);//=273+100;

  equation
  /* geometry*/
   A_a= pi*d_a*L*n;
   A_i= pi*d_i*L*n;
   A_m=(A_a-A_i)/ln(A_a/A_i);

   d_m=(d_a-d_i)/ln(d_a/d_i);

  (pi/4*d_i^2)*L*n=HU_l+HU_v;
  h_liq=HU_l/((pi/4*d_i^2)*L*n);

  /* outflow */
  vdot_l_out=1.848*h_lw*(abs(h_liq-h_w))^1.5; //outflow equation

  /* mass balance */
  der(HU_l*rho_l_out+HU_v*rho_v_out)=vdot_v_in*rho_v_in-vdot_v_out*rho_v_out-vdot_l_out*rho_l_out;

  /* energy balance */
  der(u_l*rho_l_out*HU_l+u_v*rho_v_out*HU_v)=vdot_v_in*h_v_in*rho_v_in-vdot_v_out*h_v_out*rho_v_out-vdot_l_out*h_l_out*rho_l_out-Q_kon;

  /* saturation temperature */
  T_s_out=Water.saturationTemperature(p_s);
  sat.Tsat=T_s_out;
  sat.psat=p_s;
  //der(p_s)/der(T_s_out)=2257e3/T_s_out*(HU_v-HU_l); //clausius clapeyron

  /* heat transfer */
  Q_kon=k*A_m*dT_m;

  (T_s_in-T_f_out)/(T_s_out-T_f_in)=exp(((T_s_in-T_f_out)-(T_s_out-T_f_in))/dT_m); //mean temperature difference for countercurrent flow

  k=1/((1/alpha_i)*(d_a/d_i)+(s/lambda_wall)*(d_a/d_m)+1/alpha_a);

  /* for calculation ov average coefficient of heat transfer */

  sigma_star=E*sigma;
  tau_s=sigma_star/8*rho_v_out*v^2;
  sigma=(sigma_lam^3+sigma_turb)^(1/3);
  sigma_lam=64/Re_s;
  sigma_turb=0.184/Re_s^0.2;

  E=psi/(exp(psi)-1);
  psi=-abs((Q_kon/A_i*1/h_evap)/(sigma/8*rho_v_out*v));
  L_char=(ny_water/g)^(1/3);
  Re_f=Tau/eta_water;
  Re_s=v*d_i/ny_steam;
  v=max(1e-7,vdot_v_out)/(n*(pi/4)*d_i^2);
  Tau=mdot_cond/(n*d_i*pi);
  tau_s_star=tau_s/(g*rho_v_out*(1-(rho_v_out/rho_l_out))*L);
  Ga^(1/3)=L/L_char;
  Ph=cp_v_out*(T_s_out-T_wall)/h_evap;
  Nu_turb_total*Ga^(1/3)*Ph=Re_f*Pr;
  Pr=ny_water*rho_l_out*cp_l_out/lambda_water;

  ny_steam=eta_steam/rho_v_out;
  ny_water=eta_water/rho_l_out;

  alpha_i=Nu_turb_total*lambda_water/L_char;
  alpha_a=1.95*(Q_kon/A_a)^0.72*(p_sump/1e5)^0.24;

  Nu_turb_total=((f_well*Nu_lam)^m+Nu_turb^m)^(1/m);
  Nu_turb=a*(Ph*Ga^(1/3)/Pr)^b*Pr^c*(1+e*tau_s_star^f)^d;

  Nu_lam_1=0.943*((1-rho_v_out/rho_l_out)/((Ph*Ga^(1/3))/Pr))^(1/4);
  Nu_lam_2=1.04*(((1-rho_v_out/rho_l_out)*tau_s_star)/((Ph*Ga^(1/3))/Pr))^(1/3);
  Nu_lam=(Nu_lam_1^3.5+Nu_lam_2^3.5)^(1/3.5);

  if Re_f<=1 then
    f_well=1;
  else
    f_well=Re_f^0.04;
  end if;

  if tau_s_star <= 1e-7 then
    a=2.137e-4;
    b=0.6181;
    c=0.9206;
    d=0;
    e=0;
    f=0;
    m=1.67;
  elseif tau_s_star>1e-7 and tau_s_star<=5 then
    a=2.137e-4;
    b=0.6181;
    c=0.9206;
    d=1.618;
    e=0.145;
    f=0.541;
    m=1.67;
  elseif tau_s_star>5 and tau_s_star<=10 then
    a=7.894e-3;
    b=0.2612;
    c=0.6306;
    d=1.2612;
    e=0.407;
    f=0.42;
    m=2.3;
  //elseif tau_s_star>10 and tau_s_star<=40 then
  else
    a=2.761e-2;
    b=0.1064;
    c=0.5065;
    d=1.1064;
    e=0.6469;
    f=0.473;
    m=2.75;
  end if;

  /* partial balance to gain wall temperature */
  Q_kon=alpha_i*A_i*(T_s_out-T_wall);

  initial equation
  //T_out=273+120;
  //HU_l=0.05*A*H;
  h_liq=h_w;
  //vdot_v_out=vdot_v_in;

  //assert(Q_kon> 0, "Q_kon< 0; heat is flowing in the wrong direction");

    annotation (Diagram(graphics), Icon(graphics));
  end dT_m_2;

  model alpha

    Modelica.Blocks.Sources.Ramp velocity(
      duration=100,
      height=4.9,
      offset=0.1) 
      annotation (Placement(transformation(extent={{-52,42},{-42,53}})));

  //for calculation ov average coefficient of heat transfer
  Real sigma; //coefficient of drag without condensation (1-phase)
  Real sigma_lam;
  Real sigma_turb;
  Real sigma_star; //coefficient of drag with condensation (2-phase)
  SI.ShearStress tau_s; //shear stress of steam
  SI.ShearStress tau_s_star; //normalized shear stress
  SI.Density rho_v_out=1.1;
  SI.Density rho_l_out=950;
  SI.Velocity v=velocity.y; //velocity of steam
  Real E; //         amplification factor
  Real psi;
  SI.CoefficientOfHeatTransfer alpha_i(start=1e2);// valid for laminar and turbulent flow; VDI-Wärmeatlas
  SI.CoefficientOfHeatTransfer alpha_a(start=1e4);//for bulk boiling in pressure vessels (natural convection); Stephan, Baehr, "Wärme- und Stoffübertragung"
  Real Nu_turb_total;//Nusselt number; valid for laminar and turbulent flow
  Real Nu_lam;
  Real Nu_lam_1;
  Real Nu_lam_2;
  Real Nu_turb;
  Real Re_f(start=1); //Reynolds number of the condensate film
  Real Re_s; //Reynolds number of steam
  SI.Length L=6;
  SI.Length L_char; //characteristic length
  SI.KinematicViscosity ny_water;
  SI.KinematicViscosity ny_steam;
  SI.DynamicViscosity eta_water=4.65e-4;
  SI.DynamicViscosity eta_steam=1.295e-5;
  Real Tau;
  Real Ga; //Galilei number
  Real Ph; //phase change number
  Real Pr; //Prandtl number
  Real a;
  Real b;
  Real c;
  Real d;
  Real e;
  Real f;
  Real m; //tabulated value - VDI Wärmeatlas
  Real f_well;
  parameter Real n=80; //number of tubes inside the heat exchanger
  SI.MassFlowRate mdot_cond;
  SI.ThermalConductivity lambda_water=0.66;
  constant Real g=Modelica.Constants.g_n;
    import Modelica.Constants.pi;
  SI.HeatFlowRate Q_kond=1e5;
  SI.SpecificHeatCapacity cp_l_out = 4187;
  SI.SpecificHeatCapacity cp_v_out = 2027;
  SI.SpecificEnthalpy h_evap=2.2e6;
  SI.Temperature T_wall=T_v_out-25;
  SI.Temperature T_v_out=120+273.15;
  SI.Pressure p=2e5;

  /* Geometry */
  parameter SI.Length d_a=25/1000; //outer tube diameter
  parameter SI.Length d_i=21/1000; //inner tube diameter
  parameter SI.Length s=(d_a-d_i)/2; //wall thickness of tube
  SI.Length d_m;
  SI.Area A_a; // inner surface area of the tubes
  SI.Area A_i; // outer surface area of the tubes
  SI.Area A_m; // mean value of surface area

  equation
  sigma_star=E*sigma;
  tau_s=sigma_star/8*rho_v_out*v^2;
  sigma=(sigma_lam^3+sigma_turb)^(1/3);
  sigma_lam=64/Re_s;
  sigma_turb=0.184/Re_s^0.2;

  E=psi/(exp(psi)-1);
  psi=abs((Q_kond/A_i*1/h_evap)/(sigma/8*rho_v_out*v));
  L_char=(ny_water/g)^(1/3);
  Re_f=Tau/eta_water;
  Re_s=v*d_i/ny_steam;
  Tau=mdot_cond/(n*d_i*pi);
  tau_s_star=tau_s/(g*rho_v_out*(1-(rho_v_out/rho_l_out))*L);
  Ga^(1/3)=L/L_char;
  Ph=cp_v_out*(T_v_out-T_wall)/h_evap;
  Nu_turb_total*Ga^(1/3)*Ph=Re_f*Pr;
  Pr=ny_water*rho_l_out*cp_l_out/lambda_water;

  ny_steam=eta_steam/rho_v_out;
  ny_water=eta_water/rho_l_out;

  alpha_i=Nu_turb_total*lambda_water/L_char;
  alpha_a=1.95*(Q_kond/A_a)^0.72*p^0.24;

  Nu_turb_total=((f_well*Nu_lam)^m+Nu_turb^m)^(1/m);
  Nu_turb=a*(Ph*Ga^(1/3)/Pr)^b*Pr^c*(1+e*tau_s_star^f)^d;

  Nu_lam_1=0.943*((1-rho_v_out/rho_l_out)/((Ph*Ga^(1/3))/Pr))^(1/4);
  Nu_lam_2=1.04*(((1-rho_v_out/rho_l_out)*tau_s_star)/((Ph*Ga^(1/3))/Pr))^(1/3);
  Nu_lam=(Nu_lam_1^3.5+Nu_lam_2^3.5)^(1/3.5);

  if Re_f<=1 then
    f_well=1;
  else
    f_well=Re_f^0.04;
  end if;

  if tau_s_star <= 1e-7 then
    a=2.137e-4;
    b=0.6181;
    c=0.9206;
    d=0;
    e=0;
    f=0;
    m=1.67;
  elseif tau_s_star>1e-7 and tau_s_star<=5 then
    a=2.137e-4;
    b=0.6181;
    c=0.9206;
    d=1.618;
    e=0.145;
    f=0.541;
    m=1.67;
  elseif tau_s_star>5 and tau_s_star<=10 then
    a=7.894e-3;
    b=0.2612;
    c=0.6306;
    d=1.2612;
    e=0.407;
    f=0.42;
    m=2.3;
  //elseif tau_s_star>10 and tau_s_star<=40 then
  else
    a=2.761e-2;
    b=0.1064;
    c=0.5065;
    d=1.1064;
    e=0.6469;
    f=0.473;
    m=2.75;
  end if;

  /* geometry */
   A_a= pi*d_a*L*n;
   A_i= pi*d_i*L*n;
   A_m=(A_a-A_i)/ln(A_a/A_i);

   d_m=(d_a-d_i)/ln(d_a/d_i);

  end alpha;
end HeatTransferModel;
