within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model Virial_SecondVirialCoefficient
                                     // Properties of Gases & Liquids 5th Ed. (4.13, 5.10)

  replaceable package MediumVapour = 
      ThermalSeparation.Media.C2H5OH_Water_Vap   constrainedby
    ThermalSeparation.Media.BaseMediumVapour;

    parameter Integer n=1;
    parameter Integer nS=2 "number of components";
    parameter Integer reorgVap[nS] = {1,2};

  final parameter SI.Temperature Tci[nS]= {MediumVapour.Tcrit[reorgVap[i]] for i in 1:nS};
  final parameter SI.Pressure pci[nS]= {MediumVapour.pcrit[reorgVap[i]] for i in 1:nS};
  final parameter Real omega[nS] = {MediumVapour.omega[reorgVap[i]] for i in 1:nS};
  final parameter Units.DipoleMoment my[nS] = {MediumVapour.mu[reorgVap[i]] for i in 1:nS};

  parameter SI.MolarVolume Vci_m3_mol[nS] = {MediumVapour.Vcrit[reorgVap[i]] for i in 1:nS}
    "Vektor der kritischen Volumina in m3/mol";
  input SI.Temperature T[n] "[K]";
  input Integer NoOfEq[nS];

  output Real B[n,nS,nS] "in m³/mol";

  Virial_TsonopoulosConstants TsonopoulosConstants(redeclare replaceable
      package MediumVapour = 
      MediumVapour, nS=nS, reorgVap=reorgVap, NoOfEq=NoOfEq);

protected
  Real Vci[nS] = Vci_m3_mol*1e6 "critical molar volumina in cm3/mol";
  Real Zci[nS] "compressibility factors at the critical point";
  Real kij[nS,nS];
  Real Tcik[nS,nS];
  Real Vcik[nS,nS];
  Real pcik[nS,nS];
  Real Zcik[nS,nS];
  Real omegaik[nS,nS];
  Real Tr[n,nS,nS];
  Real aik[nS,nS];
  Real bik[nS,nS];
  Real f0[n,nS,nS];
  Real f1[n,nS,nS];
  Real f2[n,nS,nS];
  Real f3[n,nS,nS];
  Real Summe[n,nS,nS];

//Tsonopoulos Parameter für 2ten Virialkoeff.
  Real ai[nS]= TsonopoulosConstants.a;
  Real bi[nS]= TsonopoulosConstants.b;

equation
for i in 1: nS loop
  Zci[i] = pci[i]*Vci_m3_mol[i]/Modelica.Constants.R/Tci[i];
end for;

for i in 1:nS loop
  for k in 1:nS loop
    if i==k then
            Tcik[i,k]=  Tci[i];
        Vcik[i,k]=  Vci_m3_mol[i];
        pcik[i,k]= pci[i];
        Zcik[i,k]= Zci[i];
        omegaik[i,k]= omega[i];
        kij[i,k]= 0;
        aik[i,k] = ai[i];
        bik[i,k] = bi[i];
    else
              Vcik[i,k]=   ((Vci_m3_mol[i]^(1/3)+Vci_m3_mol[k]^(1/3))/2)^3;
        Zcik[i,k]=   (Zci[i]+Zci[k])/2;

        Tcik[i,k]=   (Tci[i]*Tci[k])^(1/2)*(1-kij[i,k]);
        pcik[i,k]=   Zcik[i,k]*Modelica.Constants.R*Tcik[i,k]/Vcik[i,k];
        omegaik[i,k]=   (omega[i]+omega[k])/2;

        if my[i]== 0 and my[k] == 0 then
          kij[i,k]=   1-((2*(Vci[i]*Vci[k])^(1/6))  /  (Vci[i]^(1/3)+Vci[k]^(1/3))^3);
        else
          kij[i,k] = 0.154;
        end if;

        if my[i]== 0 or my[k] == 0 then
          aik[i,k] = 0;
          bik[i,k] = 0;
        else
          aik[i,k]= (ai[i]+ai[k])/2;
          bik[i,k]= (bi[i]+bi[k])/2;
        end if;
        end if;
        end for;
        end for;

for j in 1:n loop
  for i in 1:nS loop
    for k in 1:nS loop

        f0[j,i,k] = 0.1445-0.330/Tr[j,i,k]-0.1385/Tr[j,i,k]^2-0.0121/Tr[j,i,k]^3-0.000607/Tr[j,i,k]^8;
        f1[j,i,k] = 0.0637+0.331/Tr[j,i,k]^2-0.423/Tr[j,i,k]^3-0.008/Tr[j,i,k]^8;
        f2[j,i,k] = 1/Tr[j,i,k]^6;
        f3[j,i,k] = -1/Tr[j,i,k]^8;

        Summe[j,i,k] = f0[j,i,k] + omegaik[i,k]*f1[j,i,k]+aik[i,k]*f2[j,i,k]+bik[i,k]*f3[j,i,k];

      if i==k then

        Tr[j,i,k]=  T[j]/Tci[i];
        B[j,i,k] = Modelica.Constants.R*Tci[i]/pci[i]*Summe[j,i,k] "in m³/mol";
        //Bij[i,j]= ((0.083-(0.422/Tr[i,j]^(1/6)))+omega[i]*(0.139-(0.172/Tr[i,j]^(4.2))))*(Modelica.Constants.R*Tci[i])/pci[i];
      else

        Tr[j,i,k]=  T[j]/Tcik[i,k];

        B[j,i,k] = Modelica.Constants.R*Tcik[i,k]/pcik[i,k]*Summe[j,i,k]
            "in m³/mol";

        //Bij[i,j]= ((0.083-(0.422/Tr[i,j]^(1/6)))+omegaij[i,j]*(0.139-(0.172/Tr[i,j]^(4.2))))*(Modelica.Constants.R*Tcij[i,j])/pcij[i,j];
      end if;
    end for;

  end for;
end for;

// for i in 1: nS loop
//   Zci[i] = pci[i]*1e5*Vci_m3_mol[i]/Modelica.Constants.R/Tci[i];
// end for;
//
//   for i in 1:nS loop
//     for j in 1:nS loop
//       if i==j then
//         Tcij[i,j]=  Tci[i];
//         Vcij[i,j]=  Vci[i];
//         pcij[i,j]= pci[i];
//         Zcij[i,j]= Zci[i];
//         omegaij[i,j]= omega[i];
//         kij[i,j]= 0;
//         Tr[i,j]=  T/Tci[i];
//         B[m,k] = Modelica.Constants.R*10*Tcrit[k]/pcrit_bar[k]*Summe[m,k];
// //Bij[i,j]= ((0.083-(0.422/Tr[i,j]^(1/6)))+omega[i]*(0.139-(0.172/Tr[i,j]^(4.2))))*(Modelica.Constants.R*Tci[i])/pci[i];
//       else
//         Vcij[i,j]=   ((Vci[i]^(1/3)+Vci[j]^(1/3))/2)^3;
//         Zcij[i,j]=   (Zci[i]+Zci[j])/2;
//         kij[i,j]=   1-((8*(Vci[i]*Vci[j])^(1/2))  /  (Vci[i]^(1/3)+Vci[j]^(1/3))^3);
//         Tcij[i,j]=   (Tci[i]*Tci[j])^(1/2)*(1-kij[i,j]);
//         pcij[i,j]=   Zcij[i,j]*Modelica.Constants.R*Tcij[i,j]/Vcij[i,j];
//         omegaij[i,j]=   (omega[i]+omega[j])/2;
//         Tr[i,j]=  T/Tcij[i,j];
//         //Bij[i,j]= ((0.083-(0.422/Tr[i,j]^(1/6)))+omegaij[i,j]*(0.139-(0.172/Tr[i,j]^(4.2))))*(Modelica.Constants.R*Tcij[i,j])/pcij[i,j];
//       end if;
//     end for;
//   end for;
end Virial_SecondVirialCoefficient;
