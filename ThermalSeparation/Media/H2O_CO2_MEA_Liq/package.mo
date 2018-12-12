within ThermalSeparation.Media;
package H2O_CO2_MEA_Liq "H2O,CO2,MEA: CO2 separation with MEA"
extends ThermalSeparation.Media.BaseMediumLiquidReaction(has_etaSubstance={true,false,false},
    nSubstance=3,
    Tcrit={647.3,304.12,500},
    pcrit={2.2064e7,7.374e6,1e6},
    Vcrit={55.95e-6,94.07e-6,150e-6},
    omega={0.344,0.225,0.5},
    mu={1.8,0,1.5},
    MMX={0.018,0.044,0.061},
    eq_Tsonopoulos={6,1,4},
    henry=fill(false, nSubstance),
    nR=1,
    reacComp={2},
    nu={{-1,-1,-1}},
    nS_reac = 3);
   // all values except for water properties and molar masses are dummy values

constant Integer ic[nSubstance]=zeros(nSubstance);

                                                           //fill(false,nSubstance);


  redeclare replaceable model extends BaseProperties

   Real alpha; // loading in mol CO2 / mol MEA
   // Real m[nSubstance]; // molarity in mol i / kg H2O
   Real m_av; // average molarity
   Modelica.SIunits.Temperature Theta = T - 273.15;
                                      // Temperature in °C
   Real d_koeff[8];
   CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0, T=T, p=p, x=x);

  equation
      eta_comp=fill(eta, nSubstance);
   /* calculation of spec. enthalpy */

     h = calcSpecificEnthalpy.h;

   /* calculation of density (cf. Oexmann) */

      d = d_koeff[1] + d_koeff[2] * Theta + d_koeff[3] * Theta^2 + d_koeff[4] * alpha + d_koeff[5] * alpha^2 + d_koeff[6] * m_av +
          d_koeff[7] * alpha^2 + d_koeff[8] * alpha * m_av;

     //d = 986.9;

     d_koeff = {1.005e3,-5.993e-1,0,2.492e1,1.135e2,6.172e0,-3.57e-1,1.917e1};

     m_av = x[3]/(x[1]*MMX[1]);

     // m_av = sum(m[i] for i in 1:nSubstance)/nSubstance;

     //alpha = x[2] / x[3]; // c_carba = c_CO2;
     alpha = x[2]/max(1e-7,x[3]);

     MM = sum(x[i]*MMX[i] for i in 1:nSubstance);

   /* calculation of other outputs */

     u = h - p / d * MM;

     eta = 0.001002; // for water @ 20°C

     sigma = 0.0728; // for water @ 20°C

     v = MM / d;

      for i in 2:nSubstance loop
      p_sat[i] = 1e5; // dummy value
      end for;
      Modelica.Math.log(p_sat[1]) = 73.649 - 7258.2 / T - 7.3037 * Modelica.Math.log(T) + 4.1653e-6 * T^2; //water

      lambda =0.55;
      cp=4000;

  end BaseProperties;


  redeclare replaceable model extends ActivityCoefficient

  equation
  gamma = fill(1,nSubstance);

  end ActivityCoefficient;


  redeclare model extends FugacityCoefficient

  equation
  phi_sat = fill(1,nSubstance);

  end FugacityCoefficient;


  redeclare model extends HenryCoefficient

  equation
  // He=fill(1,nSubstance); // not needed or neglected

  end HenryCoefficient;


  redeclare model extends DiffusionCoefficient(diluted=true)
  /*** model which shall be used to calculate the diffusion coefficients in a binary mixture ***/

  equation
     D[1] = 2.35e-6 * Modelica.Math.exp(-2119/T); //diff. coeff H2O - CO2
     D[2] = 2.35e-6 * Modelica.Math.exp(-2119/T); //diff. coeff H2O - MEA (at the moment still same coeff as H2O - CO2)
     D[3] = 2.35e-6 * Modelica.Math.exp(-2119/T); //dummy value for diff coeff CO2 - MEA
     D_diluted[1] = D[1]; //diff. coeff at infinity dissolution of CO2 in H2O (at the moment still same coeff as H2O - CO2)
     D_diluted[2] = D[2]; //diff. coeff at infinity dissolution of MEA in H2O (at the moment still same coeff as H2O - CO2)
     D_diluted[3] = D[3]; //self- diffusion coeff. of water (at the moment still same coeff as H2O - CO2)

  end DiffusionCoefficient;


  redeclare model extends CalcSpecificEnthalpy

  Real alpha(start=0.2); // loading in mol CO2 / mol MEA
  Modelica.SIunits.SpecificHeatCapacity cp_mix;
                                  // heat capacity of solution
  Modelica.SIunits.SpecificHeatCapacity cp_tot;
                                  // introduced for better initialization
  Real cp_koeff[11];
  // Real m[nSubstance]; // molarity in mol i / kg H2O
  Real m_av; // average molarity
  Modelica.SIunits.Temperature Theta=T - 273.15;
                                     // Temperature in °C
  Modelica.SIunits.MolarMass MM;

  parameter Real omega_k=0.05
    "large value if change between constant and variable shall be steep";
  parameter Real omega_time = 500 "Wendepunkt der tanh-Funktion";
  Real omega = 0.5*tanh(omega_k*(time-omega_time))+0.5;

  equation
   /* calculation of spec. heat capacity (cf. Oexmann) */

     cp_mix / 1000 =  cp_koeff[1] + cp_koeff[2] * Theta + cp_koeff[3] * Theta^2 + cp_koeff[4] * alpha + cp_koeff[5] * alpha^2 + cp_koeff[6] * m_av +
               cp_koeff[7] * m_av^2 + cp_koeff[8] * Theta * alpha + cp_koeff[9] * Theta * m_av + cp_koeff[10] * alpha * m_av +
               cp_koeff[11] * Theta * alpha * m_av;

     cp_koeff = {4.294e0,-1.859e-3,2.575e-5,-7.819e-1,6.536e-1,-1.124e-1,4.746e-3,8.181e-4,8.267e-5,-5.364e-2,-1.909e-4};

  //      for i in 2:nSubstance loop
  //      m[i] = x[i]/(max(x[1],1e-7)*MMX[1]);
  //      end for;

     m_av = x[3]/(x[1]*MMX[1]);

     alpha = x[2] / max(x[3],1e-7);

     /* calculation of spec. enthalpy */

     cp_tot = min(4500,cp_mix);

     h = cp_tot * (T-T0) * MM;
     //h = 4100 * (T-T0) * MM;

     MM = sum(x[i]*MMX[i] for i in 1:nSubstance);

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;


  redeclare replaceable model extends ThermodynamicFactor
  equation
  Gamma =  diagonal(ones(nSubstance));
  end ThermodynamicFactor;


  redeclare model extends MolarReactionEnthalpy
  input Real alpha;
  Real h_koeff[ 4] = {-7904,-16810,26480,8295};
  equation

    h_R[1] = (h_koeff[1] + h_koeff[2] * alpha + h_koeff[3] *alpha^2 + h_koeff[4] * alpha^3);
  end MolarReactionEnthalpy;


end H2O_CO2_MEA_Liq;
