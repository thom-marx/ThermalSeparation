within ThermalSeparation.Media;
package Methylacetatsynthese_Vap "Methylacetat, Essigsäure, Methanol, Wasser"
  extends BaseMediumVapour(nSubstance=n, V = {73.78, 53.26, 31.25, 13.1},  eq_Tsonopoulos = {0, 0,0, 6},  sigma = {4.936, 4.5, 3.626, 2.641},
 epsilon_k = {469.8, 470, 481.8, 809.1}, MMX={0.0741,0.0601,0.03204,0.018});
  constant Integer n = 4 "number of substances";
  constant Integer nS=n "number of substances";
    constant Integer nX=n "number of substances";
   constant SI.Temperature Tcrit_user[n]= {506.1, 592.7, 513.1, 647.3};

       constant Integer henry[:] = {0};
     constant SI.Pressure henry_H[:] = {0,0,0,0} "Henry coefficient";
     constant Real henry_C[:]={0,0,0,0}
    "constant to calculate temperature dependency";
     constant SI.Temperature henry_T_norm=298 "norm temperature";
      constant SI.MolarMass MM_substances[n] = MMX;

   redeclare replaceable model extends BaseProperties( T0 = 298.15)

    CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0,T=T, x=x, p=p);

       Real MM_min = min(MMX);

   EvaporationEnthalpy evapEnthalpy(p=p, T=T);

     /*** saturation pressure ***/
     Real A_Ant[n] = {16.2685, 17.4068, 18.6073, 18.585};
     Real B_Ant[n] = {2665.5416, 3785.565, 3643.3113, 3984.9228};
     Real C_Ant[n] = {-53.424, -39.626, -33.4240, -39.724};
     SI.Pressure p_sat[n];

     /*** fugacity coefficients -ok ***/
     FugacityCoefficient fugacityCoefficient( T=T, p=p, x=x, v=MM/d);
     FugacityCoefficientPure fugacityCoefficientSat(nS=nS, T=T, p=p_sat, x={1,1,1,1});
     Real phi[nS] = fugacityCoefficient.phi;
     Real phi_sat[nS] = fugacityCoefficientSat.phi;
     Real phi_div[nS];// = phi./phi_sat;
     parameter Real k = 1e-4;
     Real omega = 0.5*tanh(0.05*(time-1200))+0.5;

   equation
       for i in 1:nX loop
    X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
   end for;

     /*** saturation pressure -ok ***/
   for i in 1:n loop
     p_sat[i] = 133.3*exp(A_Ant[i] - B_Ant[i]/(T+C_Ant[i]));
     phi_div[i] =  1*(1-omega)+  omega*(min(100,phi[i])/phi_sat[i]);
   end for;

   eta = 1.5e-5;
   /*** ideales Gasgesetz sehr gute Näherung im Vergleich zu den Werten aus Multiflash -ok ***/
   d = p/(Modelica.Constants.R*T)*MM;

   /*** molar mass ***/
   MM = max(MM_min,sum(x[j]*MMX[j] for j in 1:nS));

   //d=sum(c.*MMX);
   h = calcSpecificEnthalpy.h;//evaporation enthalpy included;
     u = h - Modelica.Constants.R*T;

       lambda =0.0248;
       cp=2000;

     annotation (Icon(graphics={
                           Rectangle(
               extent={{-100,100},{100,-100}},
               lineColor={0,0,255},
               fillColor={255,255,255},
               fillPattern=FillPattern.Solid)}));
   end BaseProperties;

  model Testmodell
    parameter Integer n=1;
    parameter Integer nS=4;
  parameter SI.Temperature T = 337;
  parameter SI.Pressure p = 100000;
  parameter SI.MoleFraction x[nS] = {0.025,0.025,0.75,0.2};//,0.05}};

  BaseProperties BP(T=T,p=p,x=x,c=fill(1,nS),x_star=x);
  end Testmodell;

  redeclare replaceable model extends EvaporationEnthalpy

  /*** Temperaturabhängigkeit fehlt ***/
        constant ThermalSeparation.Units.MolarEnthalpy delta_hv[nX] = {31918, 24628, 36116, 40667};
        constant Real[nX-1,4] C= {{44920,0.3684999,0,0},{40179,2.6036999,-5.0031,2.7069},{52002,0.3775,0,0}};
        constant SI.Temperature[nX] Tcrit={506.1, 592.7, 513.1, 647.3};
        SI.Temperature Tr[nX-1];
        Real tau[nX-1];
        Real Y[nX-1];

  equation
    for i in 1:nX-1 loop
      //Equation DIPPR 106 from User Manual Multiflash
      Tr[i]=T/Tcrit[i];
      tau[i]=1-Tr[i];
      Y[i] = C[i,2] + C[i,3]*Tr[i]+ C[i,4]*Tr[i]^2;
      h[i] = C[i,1]*tau[i]^Y[i];
      end for;

    h[4] = delta_hv[4]
  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));

  end EvaporationEnthalpy;

  redeclare replaceable model extends CalcSpecificEnthalpy

    /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/

    /*** Wärmekapazitäten durch Polynome 4. Ordnung austauschen ***/

  SI.MolarHeatCapacity cp_MeOAc_liq = 150;
  SI.MolarHeatCapacity cp_HOAc_liq = 650;
  SI.MolarHeatCapacity cp_MeOH_liq = 120;
  SI.MolarHeatCapacity cp_H2O_liq = 4187*MMX[4];
  SI.MolarHeatCapacity cp_liq[nS] = {cp_MeOAc_liq, cp_HOAc_liq, cp_MeOH_liq, cp_H2O_liq};

    /*** Wärmekapazitäten durch Polynome 4. Ordnung austauschen ***/
  SI.MolarHeatCapacity cp_MeOAc_vap = 930;
  SI.MolarHeatCapacity cp_HOAc_vap = 730;
  SI.MolarHeatCapacity cp_MeOH_vap = 490;
  SI.MolarHeatCapacity cp_H2O_vap = (1800)*MMX[4];
  SI.MolarHeatCapacity cp_vap[nS] = {cp_MeOAc_vap, cp_HOAc_vap, cp_MeOH_vap, cp_H2O_vap};

  /*** Standardbildungsenthalpien ***/
  ThermalSeparation.Units.MolarEnthalpy h_f[nS] = {-442790, -484089, -238572, -285830};
  ThermalSeparation.Units.MolarEnthalpy h_substance[nS];

  /*** saturation pressure ***/
    Real A_Ant[nS] = {16.2685, 17.4068, 18.6073, 18.585};
    Real B_Ant[nS] = {2665.5416, 3785.565, 3643.3113, 3984.9228};
    Real C_Ant[nS] = {-53.424, -39.626, -33.4240, -39.724};
    SI.Temperature T_sat[nX];

    /*** Verdampfungsenthalpie: Temperaturabhängigkeit fehlt ***/
    EvaporationEnthalpy evapEnthalpy(p=p,T=T);
    ThermalSeparation.Units.MolarEnthalpy delta_hv_const[nX] = {31918, 24628, 36116, 40667};
    SI.Temperature T_hv[nX]={56.9+273.15, 117.9+273.15,64.7+273.15, 100+273.15}
      "temperature for which evaporation temperature delta_hv_const is given";
      ThermalSeparation.Units.MolarEnthalpy delta_hv[nX];

      Real a1[nS] = {24.5, 27, 21.4, 23.8};
      Real a2[nS] = {335.7, 389.3, 336.1, 374};
   Real omega = 0.5*tanh(0.01*(time-200))+0.5;
  equation
    /*** saturation temperature ***/
    for i in 1:n loop
      delta_hv[i] = (1-omega)*delta_hv_const[i]+omega*evapEnthalpy.h[i];//delta_hv_const[i];//*((1-T_sat[i]/Tcrit_user[i])/(1-T_hv[i]/Tcrit_user[i]))^0.38;
      T_sat[i]=B_Ant[i]/(-Modelica.Math.log(max(1e-5,p*x[i])/133.3) + A_Ant[i])-C_Ant[i];
      //T_sat[i] = a1[i]*Modelica.Math.log(max(1e-5,x[i]))+a2[i];
    //  if i==1 then
     // h_substance[i] = h_f[i]+ cp_liq[i]*(T_sat[i]-T0) + delta_hv[i] + cp_vap[i]*(T-T_sat[i]);
     // else
       h_substance[i] = h_f[i]+ cp_liq[i]*(T-T0) + delta_hv[i];// + cp_vap[i]*(T-T_sat[i]);
       //end if;
    end for;

    h = sum(x.*h_substance);

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  redeclare model extends FugacityCoefficient
    parameter Real A = -17.374;
    parameter Real B = -7290;
    Real f0 = 1e5 "standard fugacity";
    Real K_D_star= exp(A-B/T)*p/f0;
    Real sq[nSubstance];

    Real test=2*K_D_star*(2-x[2])^2;
    Real test1=1+4*K_D_star * x[2]*(2-x[2]);
     Real omega = 0.5*tanh(0.08*(time-20))+0.5;

  equation
  //      for i in 1:nS loop
  //        sq[i]= sqrt(max(1e-8,1+4*K_D_star * x[2]*(2-x[2])));
  //      end for;
  //
  //      phi[2] = (sq[2]-1)/max(1e-8,(2*K_D_star*(2-x[2])*x[2]));
  //         for i in 1:1 loop
  //        phi[i] =  (1+4*K_D_star*(2-x[2])-sq[i])/max(1e-8,(2*K_D_star*(2-x[2])^2));
  //         end for;
  //       for i in 3:4 loop
  //        phi[i] =  (1+4*K_D_star*(2-x[2])-sq[i])/max(1e-8,(2*K_D_star*(2-x[2])^2));
  //         end for;

      for i in 1:nS loop
        sq[i]=sqrt(1+4*K_D_star * x[2]*(2-x[2]));
      end for;

      phi[2] = (1-omega)*1 + omega*(sq[2]-1)/(2*K_D_star*(2-x[2])*x[2]);
         for i in 1:1 loop
        phi[i] = (1-omega)*1 + omega*(1+4*K_D_star*(2-x[2])-sq[i])/(2*K_D_star*(2-x[2])^2);
         end for;
       for i in 3:4 loop
        phi[i] = (1-omega)*1 + omega* (1+4*K_D_star*(2-x[2])-sq[i])/(2*K_D_star*(2-x[2])^2);
         end for;
  end FugacityCoefficient;

  model FugacityCoefficientPure "fugacity coefficient for pure components"
    parameter Integer nS=4;
    output Real phi[nS];
    input SI.Temperature T;
    input SI.MoleFraction x[nS];
    input SI.Pressure p[nS];
    parameter Real A = -17.374;
    parameter Real B = -7290;
    Real f0 = 1e5 "standard fugacity";
    Real K_D_star[nS] = exp(A-B/T)*p/f0;
    Real sq[nS];
  equation
    for i in 1:nS loop
      sq[i]= sqrt(1+4*K_D_star[i] * x[2]*(2-x[2]));
    end for;

    phi[2] =  (sq[2]-1)/(2*K_D_star[2]*(2-x[2])*x[2]);
       for i in 1:1 loop
      phi[i] =  1;
       end for;
     for i in 3:4 loop
      phi[i] =  1;
       end for;
  end FugacityCoefficientPure;
  annotation (Documentation(info="<html>
  
</html>"));
end Methylacetatsynthese_Vap;
