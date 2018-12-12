within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model LeeKesler
/*Die Universelle Gaskonstante wird im folgenden in der Einheit [bar*cm3/mol*K] verwendet.
Daraus folgt der Faktor 10 für die Umrechnung aus [J/mol*K]. */
/*Der Druck und der kritische Druck wird in [bar] verwendet. Die Eingaben sind in [Pa]. 
Es folgt ein Faktor 1e5 für die Umrechnung!*/
/*Das Volumen wird in [cm3/mol] verwendet. Da die Eingaben in [m3/mol] sind, muss ein Faktor von 1e6 verwendet werden*/

extends BaseFugacityCoefficient;
//parameter Real kij_user[a]=fill(0.07,a);
//parameter Real kij[nS,nS](fixed=false);
  parameter Real kij[nS,nS]= fill(0.07,nS,nS);
  //   parameter Integer a = inter[nS-1]
 //   "number of binary interaction parameters depending on the number of substances";

/*Die Variablen y1-y6 sind Tabellenwerte, die nach dem ersten Durchlauf mit Hilfe von Trm und prm bestimmt werden müssen*/
/*** Die Tabellennummern beziehen sich auf: Reid, Prausnitz: The Properties of Gases & Liquids, 4th edition, McGraw-Hill,
die Gleichungen sind in Tabelle 5-13 (S. 145) aufgelistet. ***/
  Real y1[n];//=fill(-0.001,n) "
   //log (f/p)(0) aus Tabelle für reduzierte Temp und Druck des Gemischs!, Tabelle 5-6";
  Real y2[n];//=fill(0.0015,n) "=  {-0.014}
   // log (f/p)(1) aus Tabelle für reduzierte Temp und Druck des Gemischs!, Tabelle 5-7";
  Real y3[n];//=fill(0.01,n) "= {0.466} (H0-H/RTcm)0 aus Tabelle 5-2!";
  Real y4[n];//=fill(-0.001,n) "= {0.433} (H0-H/RTcm)1 aus Tabelle 5-3";
  Real y5[n];//=fill(0.998,n) "0.8455}
   // Z0 Zur Bestimmung des Kompressibiltätsfaktors des Gemischs, Tabelle 3-2";
  Real y6[n];//=fill(0.001,n) "={-0.0335}
   // Z1 Zur Bestimmung des Kompressibiltätsfaktors des Gemischs, Tabelle 3-3";

  Real Tcm[n] "[K]";
  Real pcm[n] "[bar]";
  Real Zcm[n];
  Real Vcm[n](start=fill(500,n)) "[cm3/mol]";
  Real Vcij[n,nS,nS] "[cm3/mol]";
  Real Tcij[n,nS,nS] "[K]";
  Real Tr[n,nS] "reduzierte Werte der Komponenten";
  Real pr[n,nS];
  Real Trm[n] "reduzierte Werte des Gemischs";
  Real prm[n];
  Real Zm[n] "nach Glg 3-3.1";
  Real omegam[n](start=zeros(n));
  Real z4[n] "ln(f/P)m des Gemischs";
  Real z1[n] "H0-H/RTcm";
  Modelica.Blocks.Tables.CombiTable2D table_residualEnthalpy0[n](each table = leeKesler_residualEnthalpy0);
  Modelica.Blocks.Tables.CombiTable2D table_residualEnthalpy1[n](each table = leeKesler_residualEnthalpy1);
  Modelica.Blocks.Tables.CombiTable2D table_Z0[n](each table = leeKesler_Z0);
  Modelica.Blocks.Tables.CombiTable2D table_Z1[n](each table = leeKesler_Z1);
  Modelica.Blocks.Tables.CombiTable2D table_fugacityPressureRatio0[n](each table = leeKesler_fugacityPressureRatio0);
  Modelica.Blocks.Tables.CombiTable2D table_fugacityPressureRatio1[n](each table = leeKesler_fugacityPressureRatio1);

protected
  Real z2[n,nS] "Erste Summe in Formel 5-8.16";
  Real z3[n,nS,nS] "Summe in Formel 5-8.17";
  Real z5[n,nS] "Zweite Summe in Formel 5-8.16";
  Real z6[n,nS] "d(omega)/(dy)";
  Real z7[n,nS,nS] "dT/dy nach Formel 5-8.17";
  Real z8[n,nS,nS] "dV/dy nach Formel 5-8.18";
  Real z9[n,nS,nS] "dP/dy nach Formel 5-8.19";

// initial equation
//    //Erstellen der kij Matrix:
//      //Erstellen einer Matrix aus den binären Interaktionsparametern
//       for i in 1:nS loop
//         for m in i:nS loop
//           if i==m then
//             //die Einträge auf der Diagonalen, werden auf irgendeinen Wert gesetzt, da eh nie benutzt
//             kij[i,m] = 0.07;
//           else
//             //die Einträge des Vektors werden an die richtigen Stellen auf der Matrix verteilt
//             kij[i,m] = kij_user[m-2+i];
//             end if;
//           end for;
//           for m in 1:i-1 loop
//             //k_12 = k_21, k_13 = k_31 etc.
//             kij[i,m] = kij[m,i];
//           end for;
//     end for;

equation
  y3= table_residualEnthalpy0.y;
    table_residualEnthalpy0.u1 = Trm;
    table_residualEnthalpy0.u2 = prm;

  y4= table_residualEnthalpy1.y;
    table_residualEnthalpy1.u1 = Trm;
    table_residualEnthalpy1.u2 = prm;

  y5= table_Z0.y;
    table_Z0.u1 = Trm;
    table_Z0.u2 = prm;

  y6= table_Z1.y;
    table_Z1.u1 = Trm;
    table_Z1.u2 = prm;

  y1 = table_fugacityPressureRatio0.y;
   table_fugacityPressureRatio0.u1 = Trm;
   table_fugacityPressureRatio0.u2 = prm;

     y2 = table_fugacityPressureRatio1.y;
   table_fugacityPressureRatio1.u1 = Trm;
   table_fugacityPressureRatio1.u2 = prm;

for m in 1:n loop
  Zm[m]=y5[m]+omegam[m]*y6[m];
  Trm[m]=T[m]/Tcm[m];
  prm[m]=(p[m]/1e5)/pcm[m];
  Zcm[m]=0.2905-0.085*omegam[m];
  Tcm[m]=(1/(Vcm[m]^(1/4)))*(sum(sum(y[m,i]*y[m,j]*  Tcij[m,i,j]  * (Vcij[m,i,j])^(1/4) for i in 1:nS) for j in 1:nS));
  Vcm[m]=sum(sum(y[m,i]*y[m,j]* Vcij[m,i,j] for i in 1:nS) for j in 1:nS);
  pcm[m]=((0.2905-0.085*omegam[m])*10*Modelica.Constants.R*Tcm[m])/Vcm[m];
  omegam[m]=sum(y[m,i]*MediumVapour.omega[i] for i in 1:nS);
  z4[m]=y1[m]*Modelica.Math.log(10)+omegam[m]*y2[m]*Modelica.Math.log(10);
  z1[m]=y3[m]+omegam[m]*y4[m];
  for i in 1:nS loop
    Tr[m,i]=T[m]/MediumVapour.Tcrit[i];
    pr[m,i]=p[m]/MediumVapour.pcrit[i];
    phi_aux[m,i]=exp(z4[m]+(z1[m]/T[m])*z2[m,i]+((Zm[m]-1)/pcm[m])*z5[m,i]-(y2[m]*Modelica.Math.log(10))*z6[m,i]);
    z2[m,i]=sum(if i==j then 0 else y[m,j]*z7[m,i,j] for j in 1:nS);
    z5[m,i]=sum(if i==j then 0 else y[m,j]*z9[m,i,j] for j in 1:nS);
    z6[m,i]=sum(if i==j then 0 else y[m,j]*(-MediumVapour.omega[i]+MediumVapour.omega[j]) for j in 1:nS);
    for j in 1:nS loop
      Vcij[m,i,j]= (1/8)*((MediumVapour.Vcrit[i]*1e6)^(1/3)+(MediumVapour.Vcrit[j]*1e6)^(1/3))^3;
      Tcij[m,i,j]= (MediumVapour.Tcrit[i]*MediumVapour.Tcrit[j])^(0.5)*(1-kij[i,j]); //eq. 5-6.29

      z3[m,i,j]=sum(y[m,l]*(Vcij[m,l,j]^(1/4)*Tcij[m,l,j]-Vcij[m,l,i]^(1/4)*Tcij[m,l,i]) for l in 1:nS);

      z7[m,i,j]=(2*z3[m,i,j]-(0.25/(Vcm[m]^(3/4)))*z8[m,i,j]*Tcm[m])/(Vcm[m]^(1/4));
      z8[m,i,j]=2*(sum(y[m,l]*Vcij[m,l,j]-Vcij[m,l,i] for l in 1:nS));
      z9[m,i,j]=pcm[m]*(-(0.085*(MediumVapour.omega[j]-MediumVapour.omega[i])/Zcm[m])+(1/Tcm[m])*z7[m,i,j]-(1/Vcm[m])*z8[m,i,j]);
      end for;
  end for;
end for;

end LeeKesler;
