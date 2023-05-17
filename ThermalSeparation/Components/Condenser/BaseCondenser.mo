within ThermalSeparation.Components.Condenser;
partial model BaseCondenser

  parameter Integer nS(min=2)=2
    "number of species which are equal in vapour and liquid phase" annotation(Dialog(group="Shell side parameters"));
  final parameter Integer nL=MediumLiquid.nSubstance -nS
    "number of additional substances which are only in liquid phase";
  final parameter Integer nV = MediumVapour.nSubstance -nS
    "number of additional substances which are only in the vapour phase";
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;

parameter Integer mapping[nS,2] = {{i,i} for i in 1:nS}
    "parameter to map the different medium vectors one to another" annotation(Dialog(group="Shell side parameters"));
parameter Boolean connectCoolant=false "set true to enable coolant connectors" annotation(Dialog(group="Shell side parameters"));
parameter Modelica.Units.SI.Pressure p_nom=0.05e5
    "linear pressure drop with nominal value" annotation(Dialog(group="Shell side parameters"));
parameter SI.VolumeFlowRate vdot_nom=0.4 "nominal volume flow rate" annotation(Dialog(group="Shell side parameters"));
replaceable package MediumVapour =
      ThermalSeparation.Media.C2H5OH_Water_Vap  constrainedby ThermalSeparation.Media.BaseMediumVapour
    "medium to be used in vapour phase"                                          annotation(choicesAllMatching=true,Dialog(group="Shell side parameters"));
    MediumVapour.BaseProperties mediumVapourIn(T=T_v_in, x=x_v_in, p=p,x_star=x_v_in, c=c_v_in);
    //MediumVapour.SpecificEnthalpy enthalpyVapIn(T=T_v_in, x=x_v_in, p_sys=p,Tstar=273+78);

  replaceable package MediumLiquid =
    ThermalSeparation.Media.C2H5OH_Water_Liq  constrainedby ThermalSeparation.Media.BaseMediumLiquid
    "medium to be used in liquid phase"                                          annotation(choicesAllMatching=true,Dialog(group="Shell side parameters"));
    MediumLiquid.BaseProperties mediumLiquid(T=T_l, x=x_l, p=p,h=h_l);
    //MediumLiquid.SpecificEnthalpy enthalpyLiq(T=T_l, x=x_l, p_sys=p,n_mol={0,0},Tstar=Tstar);

replaceable package Coolant=Modelica.Media.Water.WaterIF97_ph constrainedby Modelica.Media.Water.WaterIF97_base
                                        "medium to be used as coolant" annotation(choicesAllMatching=true,Dialog(enable= not connectCoolant, group="Tube side parameters"));

SI.HeatFlowRate Q_kon;
//SI.Temperature Tstar;

//inflow
     SI.VolumeFlowRate vdot_v_in;
     SI.MassFlowRate mdot_v_in = vdot_v_in*rho_v_in;
     SI.MoleFraction x_v_in[nSV];
     SI.Temperature T_v_in;
     SI.Concentration c_v_in[nSV];
     ThermalSeparation.Units.MolarEnthalpy h_v_in;
     SI.Density rho_v_in;
     SI.MolarMass MM_v_in;
     //SI.VolumeFlowRate vdot_w_in;
     SI.MassFlowRate mdot_coolant_in;
//outflow
    //parameter SI.MoleFraction x_l_start[nSL]={0.9,0.1};
     SI.VolumeFlowRate vdot_l_out;
     SI.MoleFraction x_l[nSL];
     SI.Concentration c_l[nSL];
     ThermalSeparation.Units.MolarEnthalpy h_l;
     ThermalSeparation.Units.MolarEnthalpy u_l = mediumLiquid.u;
     SI.Density rho_l;
     SI.MolarMass MM_l;
     SI.Temperature T_l;

//pressure at top stage
     SI.Pressure p_in;
     SI.Pressure p;
     SI.Pressure p_sat[nSL]=mediumLiquid.p_sat;

equation
  MM_v_in = mediumVapourIn.MM;
  rho_v_in = mediumVapourIn.d;
  h_v_in = mediumVapourIn.h;

  //h_l =mediumLiquid.h;
  rho_l = mediumLiquid.d;
  MM_l = mediumLiquid.MM;

end BaseCondenser;
