within ThermalSeparation.FilmModel.BaseClasses;
model MaxwellStefan_mod "modified Maxwell-Stefan equations, no film reaction"
 extends ThermalSeparation.FilmModel.BaseClasses.BaseFilmModel(
                        final EQ=false);
SI.Density rho_v[n]=propsVap.rho;
SI.MolarMass MM_v[n] = propsVap.MM;

   replaceable model ThermoEquilibrium =
      ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby
    ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium annotation (Dialog(
      tab="Propagated from Column",
      group=
          "These variables are propagated from the column model and do not have to be set by the user!",
      enable=false));

  ThermoEquilibrium thermoEquilibrium[n](x_vap_liq=x_vap_liq,
    each nS=nS, each mapping =                                                  mapping,
    redeclare replaceable package MediumVapour =  MediumVapour,
    redeclare replaceable package MediumLiquid =   MediumLiquid,
    p=p_v[1:n], T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v);

Real K[n,nS] "equilibrium constant" annotation(Dialog(__Dymola_initialDialog=true));

  /*** state variables ***/
  replaceable model StateSelection =
      ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq
                                                                                                        constrainedby
    ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq
                                                                                annotation(choicesAllMatching,Dialog(enable=enableDialog));
  StateSelection stateSelection(
     redeclare replaceable package MediumVapour =  MediumVapour,
    redeclare replaceable package MediumLiquid =   MediumLiquid,
  n=n, propsLiq=propsLiq, propsVap=propsVap, c_l=c_l, p_v= p_v[1:n]);
 MediumLiquid.ThermodynamicFactor thermoFactorL[n]( T=T_l, x=x_l);
 MediumVapour.ThermodynamicFactor thermoFactorV[n];

   MediumLiquid.CalcSpecificEnthalpy liqToFilmB[n](each T0=T_ref,   p=p_hyd[1:n], T=T_l, x=x_transfer_fromL);
 MediumLiquid.CalcSpecificEnthalpy filmToLiqB[n](each T0=T_ref,  p=p_hyd[1:n], T=T_l, x=x_transfer_toL);
  MediumVapour.CalcSpecificEnthalpy vapToFilmB[n](each T0=T_ref,                       p=p_v[1:n], T=T_v, x=x_transfer_fromV);
 MediumVapour.CalcSpecificEnthalpy filmToVapB[n](each T0=T_ref,   p=p_v[1:n], T=T_v, x=x_transfer_toV);

 //to be provided by the extending class
   SI.CoefficientOfHeatTransfer alphaLiq[n];
SI.CoefficientOfHeatTransfer alphaVap[n];
   SI.Area A_I[n]( start=400*ones(n)) "interfacial area";
   Real k_l[n,nSL,nSL];
    Real k_v[n,nSV,nSV];

  MaxwellStefanMatrix[n] maxwellStefanMatrixVap(each nS = nSV, x=x_v, dummy=k_v);
   MaxwellStefanMatrix[n] maxwellStefanMatrixLiq(each nS = nSL, x=x_l, dummy=k_l);
     Real R_l[n,nSL-1,nSL-1] = maxwellStefanMatrixLiq.matrix;
  Real R_v[n,nSV-1,nSV-1]= maxwellStefanMatrixVap.matrix;
   Real R_dash_l[n,1,nSL-1];
   Real R_dash_v[n,1,nSV-1];

protected
   Real vector_liq[ n,nSL];
   Real vector_vap[ n,nSV];

   SI.MolarFlowRate Ndot_v_transfer_cond[     n,nSV](start=fill(-0.1,n,nSV))
    "all vapour is condensing";
   SI.MolarFlowRate Ndot_v_transfer_MS[     n,nSV](start=fill(-0.1,n,nSV))
    "calculated using Maxwell-Stefan eq.";

       SI.HeatFlowRate Qdot_l_transfer[n];
 SI.HeatFlowRate Qdot_v_transfer[  n];

      SI.MolarFlowRate Ndot_l_tot[        n](start=fill(1,n));
  SI.MolarFlowRate Ndot_v_tot[        n];

    Real pot_diff[n];// difference in electrostatic potential
      final parameter Integer nSLi=if sum(abs(MediumLiquid.ic))==0 then nSL-1 else nSL-2
    "number of independent composition gradients in the liquid phase";

        /*** energy balance ***/
       SI.MolarFlowRate Ndot_fromL[  n,nSL];
   SI.MolarFlowRate Ndot_toL[  n,nSL];
  SI.MolarFlowRate Ndot_fromV[ n,nSV];
  SI.MolarFlowRate Ndot_toV[    n,nSV];
    SI.MoleFraction x_transfer_fromL[n,nSL];
    SI.MoleFraction x_transfer_toL[n,nSL];
    SI.MoleFraction x_transfer_fromV[n,nSV];
    SI.MoleFraction x_transfer_toV[n,nSV];
      ThermalSeparation.Units.MolarEnthalpy h_transfer_fromV[n]= vapToFilmB.h;
  ThermalSeparation.Units.MolarEnthalpy h_transfer_toV[n]= filmToVapB.h;
  ThermalSeparation.Units.MolarEnthalpy h_transfer_fromL[n]= liqToFilmB.h;
  ThermalSeparation.Units.MolarEnthalpy h_transfer_toL[n]= filmToLiqB.h;

equation
  /*** film is calculated steady-state ***/
  if homotopyMethod.bool_Edot_inter and homotopyMethod.useHomotopy then

    Edot_v_transfer=homotopy(actual=Edot_v_interface,simplified=homotopyMethod.Edot_v_inter);
    Edot_l_transfer=homotopy(actual=Edot_l_interface,simplified=homotopyMethod.Edot_l_inter);

else

    Edot_v_transfer=Edot_v_interface;
    Edot_l_transfer=Edot_l_interface;

end if;

 if homotopyMethod.bool_Ndot_inter and homotopyMethod.useHomotopy then

  Ndot_v_transfer=homotopy(actual=Ndot_v_interface,simplified=homotopyMethod.Ndot_v_inter);
  Ndot_l_transfer=homotopy(actual=Ndot_l_interface,simplified=homotopyMethod.Ndot_l_inter);

else

  Ndot_v_transfer=Ndot_v_interface;
  Ndot_l_transfer=Ndot_l_interface;

end if;

//   Ndot_v_interface=Ndot_v_transfer;
//   Ndot_l_interface=Ndot_l_transfer;
//   Edot_v_interface=Edot_v_transfer;
//   Edot_l_interface=Edot_l_transfer;

        for j in 1:n loop
  for i in 1:nSL loop
    vector_liq[j,i] = (-x_l[j,i] + x_l_star[j,i]);
  end for;
  for i in 1:nSV loop
    vector_vap[j,i] = (-x_v[j,i] + x_v_star[j,i]);
  end for;
  end for;

/*** modified MS-equations, summation equation at phase boundary removed ***/

    for j in 1:n loop
     if considerStartUp and startUp[j] then
       for i in 1:nSV loop
         Ndot_v_transfer[j,i] = if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];
       end for;
     else
      R_v[j,1:nSV-1,1:nSV-1]* Ndot_v_transfer[j,1:nSV-1] =Ndot_v_tot[j] * R_v[j,1:nSV-1,1:nSV-1]*   x_v[j,1:nSV-1] +  A_I[j] * propsVap[j].rho/propsVap[j].MM*thermoFactorV[j].Gamma[1:nSV-1,1:nSV-1]*vector_vap[j,1:nSV-1];
      R_dash_v[j,1,1]*Ndot_v_transfer[j,nSV] + R_dash_v[j,1,2:nSV-1]*Ndot_v_transfer[j,1:nSV-2] = Ndot_v_tot[j]*R_dash_v[j,1,1] *x_v[j,nSV] + Ndot_v_tot[j]*R_dash_v[j,1,2:nSV-1]*x_v[j,1:nSV-2] + A_I[j] * propsVap[j].rho/propsVap[j].MM*vector_vap[j,nSV];
     end if;
      Ndot_v_tot[j] = sum(Ndot_v_transfer[j,:]);
      Ndot_v_transfer_cond[j,:] = fill(0,nSV);
      Ndot_v_transfer_MS[j,:] = fill(0,nSV);
    end for;

     //Liquid side
    for j in 1:n loop
      R_l[j,1:nSL-1,1:nSL-1]*(Ndot_l_transfer[j,1:nSL-1]) =Ndot_l_tot[j] * R_l[j,1:nSL-1,1:nSL-1]*   x_l[j,1:nSL-1] + A_I[j] * (1/propsLiq[j].v* thermoFactorL[j].Gamma[1:nSL-1,1:nSL-1]*vector_liq[j,1:nSL-1] + Modelica.Constants.F/(Modelica.Constants.R*T_l[j])*c_l[j,1:nSL-1].*MediumLiquid.ic[1:nSL-1]*pot_diff[j]);
      R_dash_l[j,1,1]*Ndot_l_transfer[j,nSL] + R_dash_l[j,1,2:nSL-1]*Ndot_l_transfer[j,1:nSL-2] = Ndot_l_tot[j]*R_dash_l[j,1,1] *x_l[j,nSL] + Ndot_l_tot[j]*R_dash_l[j,1,2:nSL-1]*x_l[j,1:nSL-2] + A_I[j] * propsLiq[j].rho/propsLiq[j].MM*vector_liq[j,nSL];
      Ndot_l_tot[j] = sum(Ndot_l_transfer[j,:]);
      if sum(abs(MediumLiquid.ic))==0 then
          //no charged particles are in the solution
        pot_diff[j] =0;
      else
        sum(MediumLiquid.ic[:].*x_l_star[j,:])=0; //electroneutrality condition
        pot_diff[j]= -sum(MediumLiquid.ic[:].*(-c_l_star[j,:]+c_l[j,:]))*(Modelica.Constants.R*T_l[j])/sum(MediumLiquid.ic[:].*MediumLiquid.ic[:].*c_l[j,:])/Modelica.Constants.F;
      end if;
     end for;

      /*** energy balance at the phase boundary ***/
for j in 1:n loop
  for i in 1:nSL loop
    Ndot_fromL[j,i] = -1*min(0,Ndot_l_transfer[j,i]);
    Ndot_toL[j,i]= max(0,Ndot_l_transfer[j,i]);
    x_transfer_fromL[j,i] = Ndot_fromL[j,i]/max(1e-9,sum(Ndot_fromL[j,:]));
    x_transfer_toL[j,i] = Ndot_toL[j,i]/max(1e-9,sum(Ndot_toL[j,:]));
  end for;

  for i in 1:nSV loop
    Ndot_fromV[j,i] = -1*min(0,Ndot_v_transfer[j,i]);
      Ndot_toV[j,i] = max(0,Ndot_v_transfer[j,i]);
    x_transfer_fromV[j,i] = Ndot_fromV[j,i]/max(1e-9,sum(Ndot_fromV[j,:]));
    x_transfer_toV[j,i] = Ndot_toV[j,i]/max(1e-9,sum(Ndot_toV[j,:]));
    end for;

    Qdot_l_transfer[j] = alphaLiq[j]*A_I[j]* (T_star[j] - T_l[j]);
    Qdot_v_transfer[j] = alphaVap[j]*A_I[j]*(T_star[j] - T_v[j]);
    Edot_l_transfer[j] =   Qdot_l_transfer[j] + sum(Ndot_toL[j,:])*h_transfer_toL[j] - sum(Ndot_fromL[j,:])*h_transfer_fromL[j];
    Edot_v_transfer[j] = Qdot_v_transfer[j]- sum(Ndot_fromV[j,:])*h_transfer_fromV[j] + sum(Ndot_toV[j,:])*h_transfer_toV[j];

  end for;

  /*** phase equilibrium ***/
  for j in 1:n loop
   for i in 1:nS loop
//  K[j,i]= x_v_star_eq[j,mapping[i,1]] /x_l_star[j,mapping[i,2]];
     x_v_star[j,mapping[i,1]]= K[j,i] *x_l_star[j,mapping[i,2]];
    // K[j,i] = thermoEquilibrium[j].K[i];
   end for;
  end for;

if homotopyMethod.bool_K and homotopyMethod.useHomotopy then
    for j in 1:n loop
   for i in 1:nS loop
    K[j,i] = homotopy(actual=thermoEquilibrium[j].K[i],simplified=homotopyMethod.K[i]);
   end for;

end for;
else

  for j in 1:n loop
   for i in 1:nS loop
    K[j,i] = thermoEquilibrium[j].K[i];
   end for;

end for;
end if;

     /*** Compute R_dash ***/
          for j in 1:n loop
          for i in 1:1 loop
            for m in 1:nSV-1 loop
              if i==m then
                R_dash_v[j,i,m]= x_v[j,nSV]/k_v[j,nSV,nSV-1] + sum(x_v[j,1:nSV]./k_v[j,nSV,1:nSV]);
              else
                R_dash_v[j,i,m]= -x_v[j,nSV]*(1/k_v[j,nSV,m-1] - 1/k_v[j,nSV,nSV-1]);
              end if;
            end for;
          end for;
              end for;
           for j in 1:n loop
          for i in 1:1 loop
            for m in 1:nSL-1 loop
              if i==m then
                R_dash_l[j,i,m]= x_l[j,nSL]/k_l[j,nSL,nSL-1] + sum(x_l[j,1:nSL]./k_l[j,nSL,1:nSL]);
              else
                R_dash_l[j,i,m]= -x_l[j,nSL]*(1/k_l[j,nSL,m-1] - 1/k_l[j,nSL,nSL-1]);
              end if;
            end for;
          end for;
        end for;

  annotation (Documentation(info="<html>
<p>Instead of the &QUOT;normal&QUOT; Maxwell-Stefan approach, a modified Maxwell-Stefan approach can also be used. In this case the Maxwell-Stefan equations are supplied for nSL and nSV components rather than for nSL-1 and nSV-1 components. In addition the summation equations at the phase boundary are removed.</p>
</html>"));
end MaxwellStefan_mod;
