within ThermalSeparation.Components.HeatExchanger;
model LiquidTube2

    outer ThermalSeparation.SystemTS systemTS;
parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

replaceable package MediumLiquid = 
ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                                          annotation(choicesAllMatching);
   MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, p=p, T=T_l_in, x=x_l_in);
//     MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, p=p, T=T_l_in, x=x_l_in);

     Interfaces.LiquidPortIn liquidIn(redeclare package Medium=MediumLiquid) 
                                        annotation (Placement(transformation(
          extent={{-10,88},{10,108}},  rotation=0), iconTransformation(extent={{
            -10,88},{10,108}})));

   Interfaces.LiquidPortOut liquidOut(redeclare package Medium=MediumLiquid) 
                                          annotation (Placement(transformation(
          extent={{-10,-108},{10,-88}},
                                   rotation=0), iconTransformation(extent={{-10,
            -108},{10,-88}})));

 final parameter Integer nSL = MediumLiquid.nSubstance;

/*** Medium properties ***/
Real rho_l=mediumLiquid.d;
  SI.Concentration c_l_in[nSL](nominal = 1e4);
SI.MoleFraction x_l_in[nSL];

    SI.Temperature T_l_in;

  SI.VolumeFlowRate Vdot_l(start=1e-3, nominal = 1e3);
  SI.VolumeFlowRate Vdot_l_in;

  SI.Pressure p_out;
  SI.Pressure p_in;
  SI.Pressure p;

  /*** geometry data ***/
public
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
  final parameter SI.Height H=length_out;
//  parameter Real N_w= 20 "number of resistances in flow direction";
//  parameter SI.Length s = 0.03 "distance between the center of two tubes";
  parameter SI.Length d = 0.01 "tube diameter";
 // final parameter Real a_pL = s/d;

  parameter SI.Length height_in = -0.15;
  parameter Real zeta = 0.048;
  parameter SI.Length length_out = 0.15;
  parameter SI.Diameter d_volume = 0.025;
  Real test;

  /*** for monitoring purpose only ****/
Real delta_h;
Real delta_zeta;

equation
   //downstream
   liquidIn.p = p_in;
   liquidOut.p = p_out;
   liquidIn.c = c_l_in;
   liquidOut.c = c_l_in;
   liquidIn.Vdot = Vdot_l_in;
   liquidOut.Vdot = -Vdot_l;
   liquidIn.T = T_l_in;
   liquidOut.T = T_l_in;
   liquidIn.x = x_l_in;
   liquidOut.x = x_l_in;

  p_out=p;
Vdot_l_in = Vdot_l;
  /*** pressure loss ***/
  // Vdot_v=1;
 // p_in-p_out = 5* N_w *(eps_vap*rho_v+eps_liq*rho_l)/2*a_pL/(a_pL-1) * (Vdot_v+Vdot_l)/A;
//p_in-p_out = height_in*Modelica.Constants.g_n*rho_l + zeta*length_out/(2*d*A^2) * (Vdot_l)^2 *  rho_l;

//Vdot_l = sqrt(max(1e-8,(p_in-p_out-height_in*Modelica.Constants.g_n*rho_l)*A^2*2*d/(zeta*length_out*  rho_l)));
//zeta*(Vdot_l^2*length_out*rho_l) = max(1e-10,(p_in-p_out-height_in*Modelica.Constants.g_n*rho_l)*A^2*2*d);
zeta*(Vdot_l^2*length_out*rho_l) = (p_in-p_out-height_in*Modelica.Constants.g_n*rho_l)*A^2*2*d;
//zeta = (p_in-p_out)*A^2*2*d/(Vdot_l^2*length_out*rho_l);

 // Vdot_l = sqrt(max(1e-7,((p_in-p_out) - height_in*Modelica.Constants.g_n*rho_l)/(zeta*length_out/(2*d*A^2)*  rho_l)));
 test=((p_in-p_out) - height_in*Modelica.Constants.g_n*rho_l)/(zeta*length_out/(2*d*A^2)*  rho_l);
 delta_h=height_in*Modelica.Constants.g_n*rho_l;
 delta_zeta=zeta*length_out/(2*d*A^2) * (Vdot_l)^2 *  rho_l;

   //  Vdot_l*rho_l = Modelica.Fluid.Pipes.BaseClasses.WallFriction.Laminar.massFlowRate_dp_staticHead(dp=p_in-p_out, rho_a=rho_l, rho_b=rho_l, mu_a=mediumLiquid.eta, mu_b=mediumLiquid.eta, length=length_out,
  //    diameter=d, g_times_height_ab=Modelica.Constants.g_n*height_in, roughness=2e-5,dp_small= 2);

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={Rectangle(
          extent={{-91,40},{91,-40}},
          lineColor={0,0,0},
          pattern=LinePattern.None,
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,149,255},
          origin={1,0},
          rotation=90), Text(
          extent={{-67,28},{67,-28}},
          lineColor={0,0,0},
          textString="%name",
          origin={-5,6},
          rotation=90)}));
end LiquidTube2;
