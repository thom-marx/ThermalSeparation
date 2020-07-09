within ThermalSeparation.Components.Reboiler.BaseClasses;
partial model BaseReboiler_test2
 /*** Initialization ***/
      outer ThermalSeparation.SystemTS systemTS;

    parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));
constant Boolean delta_hv_medium = MediumVapour.delta_hv_medium;
    parameter SI.MoleFraction x_l_start[nSL]= {2e-2,1 - 2e-2 - 0.057,0.057}                                 annotation(Dialog(enable=not ini_c,  tab="Initialization"));
  parameter SI.MoleFraction x_v_start[nSV] = {1 - 0.1 - 0.001,0.0005,0.01,0.0005}
                                                                        annotation(Dialog(enable=not ini_c,  tab="Initialization"));
  parameter SI.Temperature T_vapour_start= 95+273.15 annotation(Dialog(  tab="Initialization"));
  parameter SI.Temperature T_liquid_start= 95+273.15 annotation(Dialog(  tab="Initialization"));
  final parameter SI.Temperature T_v_start =  T_vapour_start;
  final parameter SI.Temperature T_l_start =  T_liquid_start;
  parameter Real eps_liq_start = 0.071204 annotation(Dialog( tab="Initialization"));
  parameter Boolean initEQ = true annotation(Dialog( tab="Initialization"));

parameter InitOption initOption=InitOption.initEQ
    "different initialization options"
    annotation(Dialog(tab="Initialization"),Evaluate=true);
  replaceable package MediumVapour =
      ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_N2
                                                           constrainedby
    ThermalSeparation.Media.BaseMediumVapour                                                          annotation(choicesAllMatching);
  MediumVapour.BaseProperties mediumVapour(c=c_v,T0=T_ref, p=p, T=T_v, x=x_v,  x_star=x_v);

replaceable package MediumLiquid =
ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS
                                                      constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                                          annotation(choicesAllMatching);
  MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, p=p, T=T_l, x=x_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,p=p, T=T_l_in, x=x_l_in);
 parameter Integer mapping[nS,2] = {{1,2},{3,1}}
    "parameter to map the different medium vectors one to another";
parameter Boolean inertVapour[nSV] = {false,true,false,true};
parameter Boolean inertLiquid[nSL] = {false, false, true} annotation(evaluate=true);
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
  final parameter Integer nL=MediumLiquid.nSubstance-nS
    "number of additional substances which are only in liquid phase";
 final parameter Integer nV = MediumVapour.nSubstance-nS
    "number of additional substances which are only in the vapour phase";
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;

  parameter Real eps_vap_start=0.5;
  parameter Boolean fixedCirculation = false;

    /*** vapour properties ***/
  SI.Density rho_v= mediumVapour.d "density of the vapour, all components";
  SI.MolarMass MM_v( start=0.028)= mediumVapour.MM
    "molar mass of the vapour mixture ";
  ThermalSeparation.Units.MolarEnthalpy h_v = mediumVapour.h;

  SI.MolarInternalEnergy u_v(stateSelect=StateSelect.prefer)= mediumVapour.u;

  /*** liquid properties ***/
  SI.Density rho_l = mediumLiquid.d "density of the liquid, all components";
      SI.Density rho_l_in = mediumLiquidIn.d;

  SI.MolarMass MM_l(start=0.018)= mediumLiquid.MM
    "molar mass of the liquid mixture";
    SI.MolarMass MM_l_in = mediumLiquidIn.MM;
  ThermalSeparation.Units.MolarEnthalpy h_l= mediumLiquid.h;
  ThermalSeparation.Units.MolarEnthalpy h_l_in=mediumLiquidIn.h;

  SI.MolarInternalEnergy u_l(stateSelect=StateSelect.prefer) =  mediumLiquid.u;

/*** Medium properties ***/
SI.Concentration c_l[nSL](stateSelect=StateSelect.prefer);
SI.MoleFraction x_l[nSL](start=x_l_start);
  SI.Concentration c_l_in[nSL];
SI.MoleFraction x_l_in[nSL];
SI.Concentration c_v[nSV](stateSelect=StateSelect.prefer);
SI.MoleFraction x_v[nSV];
    SI.Temperature T_l_in;
  SI.Temperature T_l;
  SI.Temperature T_v(stateSelect=StateSelect.prefer);

  SI.VolumeFlowRate Vdot_v(start=1e-4);
  SI.VolumeFlowRate Vdot_l;
  SI.VolumeFlowRate Vdot_l_in(start=1e-4);

  Real eps_liq(stateSelect=StateSelect.prefer);
  Real eps_vap;

   SI.MoleFraction x_l_star[nSL](start=x_l_start);

  SI.HeatFlowRate Qdot_wall;
  SI.Pressure p_out(start=1e5);
  SI.Pressure p_in;
  SI.Pressure p;

  /*** geometry data ***/
  final parameter SI.Volume V= Modelica.Constants.pi/4*(d_HX^2 - d_tube^2*n)*length_HX;
  parameter Real zeta = 8.6;
  parameter SI.Diameter d_tube = 0.025 "outer diameter of tubes";
  parameter SI.Diameter d_HX = 0.35 "inner diameter of heat exchanger";
  parameter Integer Nw = 9 "number of rows in flow direction";
  parameter Integer n = 76 "number of tubes";
  parameter SI.Length length_HX = 0.75 "length of tubes";
  parameter SI.Length t_tube = 0.002 "wall thickness of the tubes";
    final parameter SI.Area A_HT = Modelica.Constants.pi*(2*d_tube-t_tube)/2*length_HX*n
    "mean heat transfer area";
  parameter SI.Pressure p_start=1.025e5 annotation(Dialog(tab="Initialization"));

      replaceable model PressureLoss =
      ThermalSeparation.PressureLoss.Reboiler.TubeHX                                  constrainedby
    ThermalSeparation.PressureLoss.Reboiler.BasePressureLoss annotation (
      choicesAllMatching=true);
PressureLoss pressureLoss(zeta=zeta, p_in = p_hyd[1], p_out = p_hyd[2], eps_liq = eps_liq, rho_l = rho_l, rho_v = rho_v, d_HX = d_HX, length_HX = length_HX, d_tube=d_tube, Nw=Nw);

/*** StartUp ***/

  SI.Pressure p_initial "initial pressure in column"; //has to be provided in extending class
 SI.Pressure p_hyd[2] "hydraulic pressure";

SI.Pressure wassersaeule = eps_liq*d_HX*Modelica.Constants.g_n*rho_l;

equation
/*** correlation between mole fraction x and concentration c ***/
  for i in 1:nS loop
    x_l[i] = c_l[i] * MM_l /rho_l;
    x_v[i] *rho_v= c_v[i] * MM_v;
  end for;
  for i in 1:nL loop
    x_l[i+nS] = c_l[i+nS] * MM_l /rho_l;
  end for;
  for i in 1:nV loop
    x_v[i+nS] *rho_v= c_v[i+nS] * MM_v;
  end for;

/*** StartUp ***/

  p_hyd[1] =  p;
p_hyd[2] = p_out;

  /*** mole balance ***/
   for i in 1:nSV loop
   // component balance for vapour
  V* der(eps_vap*c_v[i]) =  - Vdot_v * c_v[i];
  end for;
  for i in 1:nSL loop
    // component balance for liquid
   V* der(eps_liq*c_l[i]) = Vdot_l_in*c_l_in[i] - Vdot_l * c_l[i];
  end for;
   V* der(eps_liq*rho_l/MM_l) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l*rho_l/MM_l;

       V* der(eps_vap*rho_v/MM_v) =  -  Vdot_v*rho_v/MM_v;

  eps_liq+eps_vap=1;

    /***energy balance ***/
V * der(eps_liq*sum(c_l[:])*u_l)  = Vdot_l_in*sum(c_l_in)*h_l_in - Vdot_l*sum(c_l[:])*h_l   + Qdot_wall;
 V * der(eps_vap*sum(c_v[:])*u_v)  =  - Vdot_v*sum(c_v[:])*h_v;

//      //Liquid side
//    for i in 1:nSL loop
//      if inertLiquid[i] then
//      //  x_l[i] = x_l_star[i];
//      else
//         // x_l[i] = x_l_star[i];
//      end if;
//      end for;
      x_l[1] = x_l_star[1];
       x_l[2] = x_l_star[2];
       // x_l[3] = x_l_star[3];

           for i in 3:3 loop
     if inertLiquid[i] then
       x_l[i] = x_l_star[i];
     end if;
     end for;

  /*** thermodynamic equilibrium ***/

   for i in 1:nS loop
    // x_v[mapping[i,1]]=x_l_star[mapping[i,2]];
    // x_l[mapping[i,2]]=x_l_star[mapping[i,2]];
   end for;

  Vdot_l = Vdot_v*eps_liq/eps_vap;

 Vdot_v = pressureLoss.Vdot;

/*** system pressure ***/
p=p_in-wassersaeule;

initial equation

  for i in 1:nSL loop
     if inertLiquid[i] then
        x_l[i] = x_l_start[i];
     end if;
  end for;
  for i in 3:nSV loop
   //  x_v[i] = x_v_start[i];
  end for;

  sum(x_l[:])=1;
  sum(x_v[:])=1;

  T_v = T_v_start;
  T_l = T_l_start;

  p = p_start;

 eps_liq=eps_liq_start;

  annotation (Diagram(graphics));
end BaseReboiler_test2;
