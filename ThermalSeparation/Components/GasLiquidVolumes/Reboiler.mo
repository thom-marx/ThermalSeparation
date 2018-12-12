within ThermalSeparation.Components.GasLiquidVolumes;
model Reboiler "reboiler with connectors"
    extends Icons.Icons.Reboiler;
 extends
    ThermalSeparation.Components.GasLiquidVolumes.BaseClasses.BaseReboiler2;
      Interfaces.LiquidPortIn liquidIn(redeclare package Medium=MediumLiquid) 
                                         annotation (Placement(transformation(
           extent={{-10,-90},{10,-70}}, rotation=0), iconTransformation(extent=
            {{-10,-90},{10,-70}})));

    Interfaces.GasPortOut vapourOut(   redeclare package Medium=MediumVapour) 
                                        annotation (Placement(transformation(
           extent={{-24,68},{-4,88}},rotation=0), iconTransformation(extent={{
            -24,68},{-4,88}})));
    Interfaces.LiquidPortOut liquidOut(redeclare package Medium=MediumLiquid) 
                                           annotation (Placement(transformation(
           extent={{0,68},{20,88}}, rotation=0), iconTransformation(extent={{0,
            68},{20,88}})));
   Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
     annotation (Placement(transformation(extent={{70,-24},{90,-4}})));

equation
     //upstream
    vapourOut.p = p_out;
    vapourOut.c = c_v;
    vapourOut.Vdot = -Vdot_v;
    vapourOut.T =T_v;
    vapourOut.x =x_v;
    vapourOut.p_medium = p;
    //downstream
    liquidIn.p = p_in;
   // liquidOut.p = p_out;
    liquidIn.c = c_l_in;
    liquidOut.c = c_l;
    liquidIn.Vdot = Vdot_l_in;
    liquidOut.Vdot = -Vdot_l;
    liquidIn.T = T_l_in;
    liquidOut.T = T_l;
    liquidIn.x = x_l_in;
    liquidOut.x = x_l;

/*** heat port ***/
  Qdot_wall = heatPort.Q_flow;
  T_l=heatPort.T;

//p=p_hyd[1]-wassersaeule;
/*** StartUp ***/
  p_initial = vapourOut.p "initial pressure in column";

  annotation (Icon(graphics));
end Reboiler;
