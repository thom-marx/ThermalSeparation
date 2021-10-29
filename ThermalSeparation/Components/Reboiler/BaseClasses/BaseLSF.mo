within ThermalSeparation.Components.Reboiler.BaseClasses;
partial model BaseLSF
  extends Icons.Icons.Flash;
 /*** Initialization ***/
      outer ThermalSeparation.SystemTS systemTS;
    parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));
constant Boolean delta_hv_medium = MediumVapour.delta_hv_medium;

    parameter SI.MoleFraction x_l_start[nSL]= {2e-2,1 - 2e-2 - 0.057,0.057}                                 annotation(Dialog(enable=not ini_c and not equilibrium, tab="Initialization"));
  parameter SI.MoleFraction x_v_start[nSV] = {1 - 0.1 - 0.001,0.0005,0.01,0.0005}
                                                                        annotation(Dialog(enable=not ini_c and not equilibrium, tab="Initialization"));

  parameter SI.Temperature T_vapour_start= 95+273.15 annotation(Dialog(enable = not equilibrium, tab="Initialization"));
  parameter SI.Temperature T_liquid_start= 95+273.15 annotation(Dialog(enable = not equilibrium, tab="Initialization"));
  final parameter SI.Temperature T_v_start =  T_vapour_start;
  final parameter SI.Temperature T_l_start =  T_liquid_start;
  parameter Real eps_liq_start = 0.071204 annotation(Dialog(enable = not equilibrium, tab="Initialization"));
  parameter Boolean initEQ = true annotation(Dialog(enable = not equilibrium, tab="Initialization"));

  replaceable package MediumVapour =
      ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_N2
                                                           constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                                                                      annotation(choicesAllMatching);
  MediumVapour.BaseProperties mediumVapour(c=c_v, T0=T_ref, p=p, T=T_v, x=x_v, x_star=x_v);
  MediumVapour.CalcSpecificEnthalpy vapToFilmB(T0=T_ref, p=p, T=T_v, x=x_transfer_fromV);
  MediumVapour.CalcSpecificEnthalpy filmToVapB(T0=T_ref, p=p, T=T_v, x=x_transfer_toV);

replaceable package MediumLiquid =
ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS
                                                      constrainedby ThermalSeparation.Media.BaseMediumLiquid
                                                                                                      annotation(choicesAllMatching);
  MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, p=p, T=T_l, x=x_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,p=p, T=T_l_in, x=x_l_in);
    MediumLiquid.BaseProperties mediumLiquidStar(T0=T_ref, p=p, T=T_star, x=x_l_star);
    MediumLiquid.CalcSpecificEnthalpy liqToFilmB( T0=T_ref,p=p, T=T_l, x=x_transfer_fromL);
MediumLiquid.CalcSpecificEnthalpy filmToLiqB(T0=T_ref, p=p, T=T_l, x=x_transfer_toL);
 MediumLiquid.ActivityCoefficient activityCoeff(redeclare model ActivityCoeff =
        ApplicationsThermalSeparation.Media.Correlations.ActivityCoefficient.DesorberGamma,
                                                                                  T=T_star,x_l=x_l_star);
   MediumLiquid.FugacityCoefficient fugacityCoeffSat(T=T_star, p=p, p_sat=p_sat);

 parameter Integer mapping[nS,2] = {{1,2},{3,1}}
    "parameter to map the different medium vectors one to another";
parameter Boolean inertVapour[nSV] = {false,true,false,true};
parameter Boolean inertLiquid[nSL] = {false, false, true};
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
  parameter Integer nL=1
    "number of additional substances which are only in liquid phase";
  parameter Integer nV = 2
    "number of additional substances which are only in the vapour phase";
  final parameter Integer nSL = nS + nL;
  final parameter Integer nSV = nS + nV;

  parameter Real eps_vap_start=0.5;
  parameter Boolean fixedCirculation = false;

  SI.Pressure wassersaeule = eps_liq*h_Flash*Modelica.Constants.g_n*rho_l;

// replaceable model PressureLoss =
//       ThermalSeparation.PressureLoss.Reboiler.TubeHX                                  constrainedby
//     ThermalSeparation.PressureLoss.Reboiler.BasePressureLoss annotation (
//       choicesAllMatching=true);
// PressureLoss pressureLoss(zeta=zeta, p_in = p_hyd[1], p_out = p_hyd[2], eps_liq = eps_liq, rho_l = rho_l, rho_v = rho_v, d_HX = d_HX, length_HX = length_HX, d_tube=d_tube, Nw=Nw);

//connectors
      Interfaces.LiquidPortIn liquidIn(redeclare package Medium=MediumLiquid)
                                         annotation (Placement(transformation(
           extent={{-10,-90},{10,-70}}, rotation=0), iconTransformation(extent=
            {{-10,-90},{10,-70}})));

    Interfaces.GasPortOut vapourOut( redeclare package Medium=MediumVapour)
                                        annotation (Placement(transformation(
           extent={{-24,70},{-4,90}},rotation=0), iconTransformation(extent={{
            -24,70},{-4,90}})));
    Interfaces.LiquidPortOut liquidOut(redeclare package Medium=MediumLiquid)
                                           annotation (Placement(transformation(
           extent={{-2,70},{18,90}},rotation=0), iconTransformation(extent={{-2,
            70},{18,90}})));
   Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
     annotation (Placement(transformation(extent={{70,-24},{90,-4}})));

    /*** vapour properties ***/
  SI.Density rho_v= mediumVapour.d "density of the vapour, all components";
  SI.MolarMass MM_v( start=0.028, stateSelect=StateSelect.default)= mediumVapour.MM
    "molar mass of the vapour mixture ";
  ThermalSeparation.Units.MolarEnthalpy h_v = mediumVapour.h;
    ThermalSeparation.Units.MolarEnthalpy h_transfer_fromV= vapToFilmB.h;
  ThermalSeparation.Units.MolarEnthalpy h_transfer_toV= filmToVapB.h;
  SI.MolarInternalEnergy u_v(stateSelect=StateSelect.prefer)= mediumVapour.u;

  /*** liquid properties ***/
  SI.Density rho_l = mediumLiquid.d "density of the liquid, all components";
      SI.Density rho_l_in = mediumLiquidIn.d;
        SI.Density rho_l_star= mediumLiquidStar.d;
  SI.MolarMass MM_l(start=0.018)= mediumLiquid.MM
    "molar mass of the liquid mixture";
    SI.MolarMass MM_l_in = mediumLiquidIn.MM;
      SI.MolarMass MM_l_star= mediumLiquidStar.MM;
  ThermalSeparation.Units.MolarEnthalpy h_l= mediumLiquid.h;
  ThermalSeparation.Units.MolarEnthalpy h_l_in=mediumLiquidIn.h;
    ThermalSeparation.Units.MolarEnthalpy h_transfer_fromL= liqToFilmB.h;
  ThermalSeparation.Units.MolarEnthalpy h_transfer_toL= filmToVapB.h;
  SI.MolarInternalEnergy u_l(stateSelect=StateSelect.prefer) =  mediumLiquid.u;

/*** Medium properties ***/
SI.Concentration c_l[nSL](stateSelect=StateSelect.prefer) annotation(Dialog(group="Initialization",showStartAttribute=true));
SI.MoleFraction x_l[nSL];
  SI.Concentration c_l_in[nSL];
SI.MoleFraction x_l_in[nSL];
SI.Concentration c_v[nSV](stateSelect=StateSelect.prefer) annotation(Dialog(group="Initialization",showStartAttribute=true));
SI.MoleFraction x_v[nSV];
  SI.Concentration c_l_star[nSL];
    SI.Temperature T_l_in;
  SI.Temperature T_l;
  SI.Temperature T_v(stateSelect=StateSelect.default);

  SI.VolumeFlowRate Vdot_v(start=1e-4);
  SI.VolumeFlowRate Vdot_l;
  SI.VolumeFlowRate Vdot_l_in(start=1e-4);

  Real eps_liq(stateSelect=StateSelect.prefer);
  Real eps_vap;

    SI.MolarFlowRate Ndot_v_transfer[      nSV](start=fill(-0.1,nSV));
  SI.MolarFlowRate Ndot_l_transfer[        nSL](start=fill(0.1,nSL));

//Variablen zur Enthalpiestromberechnung über die Phasengrenze
  SI.MolarFlowRate Ndot_fromL[nSL];
   SI.MolarFlowRate Ndot_fromV[nSV];
    SI.MoleFraction x_transfer_fromL[nSL];
    SI.MoleFraction x_transfer_toL[nSL];
    SI.MoleFraction x_transfer_fromV[nSV];
    SI.MoleFraction x_transfer_toV[nSV];
    SI.HeatFlowRate Edot_l_transfer;
  SI.HeatFlowRate Edot_v_transfer;
  SI.HeatFlowRate Qdot_l_transfer;
    SI.HeatFlowRate Qdot_v_transfer;
  SI.Temperature T_star(start=T_l);
SI.SpecificEnthalpy r_water = -2462.5 * T_star +3177.8e3 "evaporation enthalpy";

   SI.MoleFraction x_l_star[nSL](start=x_l_start);
      SI.MoleFraction x_v_star[nSV](start=x_v_start);

  SI.HeatFlowRate Qdot_wall;
  SI.Pressure p_out(start=1e5);
  SI.Pressure p_in;
  SI.Pressure p;
  SI.Pressure p_sat[nSL] = {mediumLiquid.p_sat[i] for i in 1:nSL};

Boolean bool_eps;

 replaceable model ThermoEquilibrium =
       ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
                                                            annotation(choicesAllMatching=true);
   ThermoEquilibrium thermoEquilibrium(nS=nS,
      mapping =                                                                                mapping, redeclare replaceable package
                          MediumVapour =
         MediumVapour,                                                                                                    redeclare replaceable package
                          MediumLiquid =
     MediumLiquid, p=p, T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v/rho_v);

        Real K[nS] "equilibrium constant";

  /*** geometry data ***/
  final parameter SI.Volume V= Modelica.Constants.pi/4*d_Flash^2*h_Flash;
//parameter Real zeta =  8.6;
  parameter SI.Diameter d_Flash=5.7;
  parameter SI.Height h_Flash = 17;
  parameter SI.Diameter d_tube = 0.025 "outer diameter of tubes";
  parameter SI.Diameter d_HX = 0.35 "inner diameter of heat exchanger";
  parameter Integer Nw = 9 "number of rows in flow direction";
  parameter SI.Length length_HX = 0.75 "length of tubes";
  parameter SI.Pressure p_start=1.025e5 annotation(Dialog(tab="Initialization"));

     //Model for reaction
  replaceable model Reaction =
      ApplicationsThermalSeparation.Reaction.CO2_Siemens            constrainedby ThermalSeparation.Reaction.BaseReaction
                                                                                                                        annotation(choicesAllMatching=true);
  Reaction reaction(propsLiq=mediumLiquid.properties,
   final n=1,final nS= nSL, c=c_l, V=V*eps_liq, Ndot_l_transfer=Ndot_l_transfer,  gamma=activityCoeff.gamma,
  redeclare package MediumLiquid =  MediumLiquid);

  /*** for monitoring purpose only ****/
  Real sum_x = sum(x_l);
  Real sum_y=sum(x_v);
  SI.VolumeFlowRate Vdot_ges = Vdot_v + Vdot_l;
  SI.MassFlowRate mdot_v= Vdot_v*rho_v;
  SI.MassFlowRate mdot_l = Vdot_l *rho_l;
  SI.MassFlowRate mdot_out=mdot_v+mdot_l;
  SI.MassFlowRate mdot_l_in = Vdot_l_in*mediumLiquidIn.d;
  SI.Density rho_mix=eps_vap*rho_v + eps_liq*rho_l;
  SI.Volume V_liq = V*eps_liq;
  SI.MolarFlowRate Ndot_v= Vdot_v * rho_v/MM_v "total molar flow rate vapour";
  SI.MolarFlowRate Ndot_l= Vdot_l * rho_l/MM_l "total molar flow rate vapour";
  SI.MolarFlowRate Ndot_l_in= Vdot_l_in * rho_l_in/MM_l_in
    "total molar flow rate vapour";

      replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry                constrainedby ThermalSeparation.Geometry.StructuredPackedColumn.Geometry annotation (
      choicesAllMatching);
Geometry geometry(n=1);

SI.HeatFlowRate Qdot_evap = if delta_hv_medium then 0 else Ndot_v_transfer[1] *0.018*r_water;

  Modelica.Blocks.Interfaces.RealOutput level_rel=eps_liq
    annotation (Placement(transformation(extent={{80,58},{122,100}})));
  Modelica.Blocks.Interfaces.RealOutput p_v=p
    annotation (Placement(transformation(extent={{80,14},{122,56}})));
/*** StartUp ***/
  SI.Pressure p_hyd[2] "hydraulic pressure";

equation
       //upstream
    vapourOut.p = p_out;
    vapourOut.c = c_v;
    vapourOut.Vdot = -Vdot_v;
    vapourOut.T =T_v;
    vapourOut.x =x_v;
    vapourOut.p_medium=p;
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

/*** system pressure ***/
//p_in=p;
p=p_in-wassersaeule;

  for i in 1:nSL loop
    //c_l_star[j,i] = x_l_star[j,i] / MM_l_star[j]*rho_l_star[j];
    c_l_star[i] = (x_l_star[i]+x_l[i])/2 / MM_l_star*rho_l_star;
  end for;

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

  /*** to do: startUp ***/
  p_hyd[1]=p;
  p_hyd[2]=p_out;

  /*** mole balance ***/
   for i in 1:nSV loop
   // component balance for vapour
  V* der(eps_vap*c_v[i]) =  - Vdot_v * c_v[i] + Ndot_v_transfer[i]; //+ Vdot_v_in*c_v_in[i];
  end for;
  for i in 1:nSL loop
    // component balance for liquid
   V* der(eps_liq*c_l[i]) = Vdot_l_in*c_l_in[i] - Vdot_l * c_l[i] + Ndot_l_transfer[i] + reaction.Ndot[i];
  end for;
  // total mole balance for liquid and vapour
  //reaction.Ndot muß nicht berücksichtigt werden, da die Reaktion ja in der Flüssigphase stattfindet und nur "Flüssigkeit zu Flüssigkeit" wird
   if eps_liq<1e-5 and Vdot_l_in<1e-8 then
     bool_eps=true;
     der(eps_liq)=0;
 else
   bool_eps=false;
    V* der(eps_liq*rho_l/MM_l) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l*rho_l/MM_l + sum(Ndot_l_transfer[:]) + sum(reaction.Ndot[:]);
    end if;
       V* der(eps_vap*rho_v/MM_v) =  -  Vdot_v*rho_v/MM_v + sum(Ndot_v_transfer[:]);// + Vdot_v_in*mediumVapourIn.d/mediumVapourIn.MM;

  eps_liq+eps_vap=1;

    /***energy balance ***/
V * der(eps_liq*sum(c_l[:])*u_l)  = Vdot_l_in*sum(c_l_in)*h_l_in - Vdot_l*sum(c_l[:])*h_l  + Edot_l_transfer  +reaction.deltaH_R + Qdot_wall -Qdot_evap;
 V * der(eps_vap*sum(c_v[:])*u_v)  =  - Vdot_v*sum(c_v[:])*h_v +Edot_v_transfer;// + Qdot_wall;// + Vdot_v_in*sum(c_v_in[:])*mediumVapourIn.h;

 //zur Berechnung der Enthalpieströme über die Phasengrenze
  for i in 1:nSL loop
    Ndot_fromL[i] = -1*min(0,Ndot_l_transfer[i]);
    x_transfer_fromL[i] = Ndot_fromL[i]/max(1e-5,sum(Ndot_fromL[:]));
       if inertLiquid[i] then
       x_transfer_toL[i] = 0;
     end if;
  end for;
     for i in 1:nS loop
       x_transfer_toL[mapping[i,2]]=  x_transfer_fromV[mapping[i,1]];
     end for;
  for i in 1:nSV loop
    Ndot_fromV[i] = -1*min(0,Ndot_v_transfer[i]);
    x_transfer_fromV[i] = Ndot_fromV[i]/max(1e-5,sum(Ndot_fromV[:]));
         if inertVapour[i] then
       x_transfer_toV[i] = 0;
     else
     end if;
    end for;
     for i in 1:nS loop
       x_transfer_toV[mapping[i,1]]=  x_transfer_fromL[mapping[i,2]];
     end for;

/*** mass transport ***/
   for i in 1:nSV loop
      if inertVapour[i] then
       Ndot_v_transfer[i] = 0;
      end if;
      end for;
      for i in 1:nSV-1 loop
    Ndot_v_transfer[i] =  - 1e6*(x_v[i] - x_v_star[i]);
    end for;
     //Liquid side
   for i in 1:nSL loop
     if inertLiquid[i] then
       Ndot_l_transfer[i] = 0;
     end if;
     end for;
     for i in 1:nSL-1 loop
     Ndot_l_transfer[i] =  - 1e6* (x_l[i] - x_l_star[i]);
   end for;

/*** summation equations at the phase boundary ***/
    sum(x_v_star)=1;
    sum(x_l_star)=1;

/*** PHASE BOUNDARY ***/
/*** mole balance at phase boundary (steady-state, including reaction term in liquid film ***/
for i in 1:nS loop
  Ndot_v_transfer[mapping[i,1]] + Ndot_l_transfer[mapping[i,2]]  = 0;
end for;

/*** energy balance at phase boundary ***/
    -Edot_v_transfer - Edot_l_transfer  =0;//-  Qdot_evap =0;

    Qdot_l_transfer =  1e8* (T_star - T_l);
    Qdot_v_transfer= 1e8*(T_star - T_v);

      Edot_l_transfer =   Qdot_l_transfer + sum(Ndot_fromV[:])*h_transfer_toL - sum(Ndot_fromL[:])*h_transfer_fromL;
      Edot_v_transfer = Qdot_v_transfer - sum(Ndot_fromV[:])*h_transfer_fromV + sum(Ndot_fromL[:])*h_transfer_toV;

  /*** thermodynamic equilibrium ***/

 for i in 1:nS loop
       x_v_star[mapping[i,1]]= K[i] *x_l_star[mapping[i,2]];
    K[i] = thermoEquilibrium.K[i];
 end for;
// entfernt: stattdessen: Vdot_l wird über Boundary vorgegeben
//  Vdot_l = Vdot_v*eps_liq/eps_vap;

/*** Der Gleichung für Druckverlust bei Umströmung von Rohrbündeln nachempfunden (VDI-Wärmeatlas), zeta wird allerdings als Parameter vorgegeben
     Sowieso alles nur mig gut, weil der Volumenstrom sich ja signifikant über die HX-Höhe ändert ***/
//Vdot_v = pressureLoss.Vdot;// sqrt(max(1e-8,(p - p_out - eps_liq*d_HX*Modelica.Constants.g_n*rho_l)*2*length_HX^2*(d_HX-10*d_tube)^2/(zeta*Nw*(eps_liq*rho_l + eps_vap*rho_v))));

initial equation
   if initEQ then

          for i in 1:nSL loop
            if inertLiquid[i] then
              x_l[i] = x_l_start[i];
            end if;
            end for;

       for i in 3:nSV loop
        x_v[i] = x_v_start[i];
       end for;
               sum(x_l[:])=1;
               sum(x_v[:])=1;

     T_v = T_v_start;
     T_l = T_l_start;
    // der(T_v)=0;
   //  der(T_l)=0;
  // der(p)=0;
   p = p_start;

Ndot_v_transfer=zeros(nSV);

 eps_liq=eps_liq_start;

    else

       for i in 1:nSL loop
       c_l[i]= x_l_start[i] / MM_l *rho_l;
       end for;
       for i in 1:nSV loop
        x_v[i] = x_v_start[i];

       end for;
     T_v = T_v_start;
     T_l = T_l_start;
 eps_liq=eps_liq_start;
  //  p_v[1:n] = p_v_start;
end if;

  annotation (Diagram(graphics), Icon(graphics));
end BaseLSF;
