within ThermalSeparation.Components.Pumps;
package OldPumps

  model idealPump
   extends Icons.Icons.Pumps;
      outer ThermalSeparation.SystemTS systemTS;
  parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

  replaceable package MediumLiquid =
  ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby Media.BaseMediumLiquid             annotation(choicesAllMatching);
    MediumLiquid.BaseProperties mediumLiquid(T0=T_ref);
      MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref);

       Interfaces.LiquidPortIn liquidIn(redeclare package Medium=MediumLiquid)
                                          annotation (Placement(transformation(
            extent={{0,-94},{20,-74}},   rotation=0), iconTransformation(extent={{-20,
              -114},{20,-74}})));

     Interfaces.LiquidPortOut liquidOut(redeclare package Medium=MediumLiquid)
                                            annotation (Placement(transformation(
            extent={{0,102},{20,122}},
                                     rotation=0), iconTransformation(extent={{-20,82},
              {20,122}})));

  parameter SI.Time T = 10 "Delay time for outlet volume flow rate V_flow_Out";
  parameter SI.TemperatureDifference dT = 0.5 "Temperature increase T_Out - T_In";
  parameter SI.VolumeFlowRate V_flow_start=0.000305
      "Fixed initial value for outlet volume flow rate V_flow_Out";

    final parameter Integer nSL = MediumLiquid.nSubstance;
    parameter SI.Concentration c_l_start[nSL]={5,50000};

  /*** Medium properties ***/
  MediumLiquid.Density rho_l;
  SI.Concentration c_l[nSL](start = c_l_start, nominal = 1e4);
  SI.MoleFraction x_l[nSL];
    ThermalSeparation.Units.MolarEnthalpy h_l;
    ThermalSeparation.Units.MolarEnthalpy u_l;
    SI.MolarMass MM_l;
    SI.Concentration c_l_in[nSL](nominal = 1e4);
  SI.MoleFraction x_l_in[nSL];
    ThermalSeparation.Units.MolarEnthalpy h_l_in;
    MediumLiquid.Temperature T_l_in;
    MediumLiquid.Temperature T_l(start=350);
    SI.VolumeFlowRate Vdot_l(start=1e-3, nominal = 1e-3);
    SI.VolumeFlowRate Vdot_l_in(start=1e-3, nominal = 1e-3);
    MediumLiquid.AbsolutePressure p_out;
    MediumLiquid.AbsolutePressure p_in;
    MediumLiquid.AbsolutePressure p(start=2e5);
    SI.VolumeFlowRate V_flow_delayed(stateSelect=StateSelect.default);

  /*** monitoring***/
  SI.MassFlowRate mdot=Vdot_l*rho_l;

  initial equation
  V_flow_delayed = V_flow_start;

  equation

    der(V_flow_delayed) = (Vdot_l-V_flow_delayed)/T;
     //downstream
     liquidIn.p = p_in;
     liquidOut.p = p_out;
     liquidIn.c = c_l_in;
     liquidOut.c = c_l;
     liquidIn.Vdot = Vdot_l_in;
     liquidOut.Vdot = -Vdot_l;//-V_flow_delayed;
     liquidIn.T = T_l_in;
     liquidOut.T = T_l;
     liquidIn.x = x_l_in;
     liquidOut.x = x_l;

    mediumLiquidIn.p=p;
    mediumLiquidIn.T=T_l_in;
    mediumLiquidIn.x=x_l_in;
    h_l_in =mediumLiquidIn.h;

    mediumLiquid.p=p;
    mediumLiquid.T=T_l;
    mediumLiquid.x=x_l;
    h_l =mediumLiquid.h;
    u_l =mediumLiquid.u;
     rho_l = mediumLiquid.d;
     MM_l = mediumLiquid.MM;

    /*** correlation between mole fraction x and concentration c ***/
    for i in 1:nSL loop
      x_l[i] = c_l[i] * MM_l /rho_l;
    end for;

      /***energy balance ***/
   //0  = Vdot_l_in * sum(c_l_in[:])*h_l_in - Vdot_l*sum(c_l[:])*h_l;
    T_l = T_l_in+dT;

      /*** mass balance ***/
     // for i in 1:nSL loop
            // A*H* der(c_l[i]) = Vdot_l_in*c_l_in[i]  - Vdot_l * c_l[i];
            0= Vdot_l_in*sum(c_l_in)- Vdot_l * sum(c_l);
  //Vdot_l=1e-4;
     // end for;

     x_l = x_l_in;
  //     sum(x_l)=1;

    p_out=p;

  //Vdot_l_in*k = p_out-p_in;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                     graphics));
  end idealPump;
end OldPumps;
