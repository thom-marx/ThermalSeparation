within ThermalSeparation.FilmModel.BaseClasses;
model FilmDiscretization "liquid film discretisized"
  extends BaseFilmModel( final EQ=false);
  parameter Integer m(min=1)=1 "number of discrete film elements";
    parameter Integer e=1 "assymetrical grid distribution for e > 1";
  SI.Length t[n] "film thickness";
  SI.Length t_z[n,m] "thickness of one film element";
  SI.Length t_temp[n] "temperature film thickness";
  SI.Length t_temp_z[n,m] "temperature film thickness of one film element";

/*** mole balances ***/
  /*** mole balance vapour film ***/
  MediumVapour.ThermodynamicFactor thermoFactorV[n];
  SI.MoleFraction vector_vap[ n,nSV];
  SI.MolarFlowRate Ndot_v_tot[n];
  /*** mole balance liquid film ***/
  MediumLiquid.ThermodynamicFactor thermoFactorL[n,m]( T=T_z[:,1:m], x=x_l_z[:,1:m,:]);
  SI.MolarFlowRate Ndot_z[n,m+1,nSL] "inter-film molar flow rates";
  SI.MoleFraction delta_x[n,m,nSL] " inter-film mole differences";
  SI.MoleFraction x_l_z[n,m+1,nSL] "inter-film mole fractions";
  Real pot_diff_z[n,m] "difference in electrostatic potential";

  final parameter Integer nSLi= if sum(abs(MediumLiquid.ic))==0 then nSL-1 else nSL-2
    "number of independent composition gradients in the liquid phase";

    //to be provided by extending class
  SI.Area A_I[n]( start=400*ones(n)) "interfacial area";
    SI.CoefficientOfHeatTransfer alphaVap[n];
    Real k_v[n,nSV,nSV];
   SI.DiffusionCoefficient D_liq_matrix[n,m,nSL,nSL] = diffCoeffLiqMatrix.D_matrix;
MediumLiquid.DiffusionCoefficient[n,m] diffCoeffLiqMatrix(
    T=T_z[:,1:m],
    p=p_film[:,1:m],
    x=x_l_z[:,1:m,:]);

   MaxwellStefanMatrix[n,m] maxwellStefanMatrixLiq(each nS = nSL, x=x_l_z[:,1:m,:], dummy=D_liq_matrix);
  MaxwellStefanMatrix[n] maxwellStefanMatrixVap(each nS = nSV, x=x_v, dummy=k_v);
  Real R_v[n,nSV-1,nSV-1]= maxwellStefanMatrixVap.matrix;
     Real B[n,m,nSL-1,nSL-1]= maxwellStefanMatrixLiq.matrix
    "Maxwell-Stefan diffusion matrix, see Taylor & Krishna, p. 25";

    replaceable model ThermoEquilibrium =
      ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
                                                            annotation (Dialog(
      tab="Propagated from Column",
      group=
          "These variables are propagated from the column model and do not have to be set by the user!",
      enable=false));

  ThermoEquilibrium thermoEquilibrium[n](x_vap_liq=x_vap_liq,
    each nS=nS,  each mapping =                                                  mapping,
    redeclare replaceable package MediumVapour =  MediumVapour,
    redeclare replaceable package MediumLiquid =  MediumLiquid,
    p=p_v[1:n], T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v);

Real K[n,nS] "thermodynamic equilibrium constant" annotation(Dialog(__Dymola_initialDialog=true));
SI.Density rho_v[n]=propsVap.rho;
SI.MolarMass MM_v[n] = propsVap.MM;

/*** energy balance ***/
  /*** energy balance vapour film ***/
  MediumVapour.CalcSpecificEnthalpy vapToFilmB[n](each T0=T_ref,  p=p_v[1:n], T=T_v, x=x_transfer_fromV);
  MediumVapour.CalcSpecificEnthalpy filmToVapB[n](each T0=T_ref,   p=p_v[1:n], T=T_v, x=x_transfer_toV);
  SI.HeatFlowRate Qdot_v_transfer[n];
  SI.MolarFlowRate Ndot_fromV[ n,nSV];
  SI.MolarFlowRate Ndot_toV[    n,nSV];
  SI.MoleFraction x_transfer_fromV[n,nSV];
  SI.MoleFraction x_transfer_toV[n,nSV];
  ThermalSeparation.Units.MolarEnthalpy h_transfer_fromV[n]= vapToFilmB.h;
  ThermalSeparation.Units.MolarEnthalpy h_transfer_toV[n]= filmToVapB.h;
  /*** energy balance liquid film ***/
  Real Edot_z[n,m+1];
  SI.Temperature T_z[n,m+1];
  SI.ThermalConductivity lambda[n] =  propsLiq.lambda;
  SI.MolarFlowRate Ndot_leftToRight[ n,m,nSL];
  SI.MolarFlowRate Ndot_rightToLeft[ n,m,nSL];
  SI.MoleFraction x_transfer_leftToRight[n,m,nSL];
  SI.MoleFraction x_transfer_rightToLeft[n,m,nSL];
  ThermalSeparation.Units.MolarEnthalpy h_leftToRight[n,m]= leftToRight.h;
  ThermalSeparation.Units.MolarEnthalpy h_rightToLeft[n,m]= rightToLeft.h;
   ThermalSeparation.Units.MolarEnthalpy h_l_z[n,m+1]= specificEnthalpy.h;

  MediumLiquid.CalcSpecificEnthalpy specificEnthalpy[n,m+1](each T0=T_ref,   p=p_film, T=T_z, x=x_l_z);
  MediumLiquid.CalcSpecificEnthalpy leftToRight[n,m](each T0=T_ref,   p=p_film[:,1:m], T=T_z[:,1:m], x=x_transfer_leftToRight);
  MediumLiquid.CalcSpecificEnthalpy rightToLeft[n,m](each T0=T_ref,  p=p_film[:,1:m], T=T_z[:,1:m], x=x_transfer_rightToLeft);

  /*** reaction ***/
  MediumLiquid.BaseProperties mediumLiquid[n,m+1](each T0=T_ref,  p=p_film, T=T_z, x= x_l_z);

  Reaction reaction[n,m](redeclare record Geometry=BaseGeometry, redeclare each replaceable package
                               MediumLiquid =
        MediumLiquid,
        propsLiq= mediumLiquid[:,1:m].properties,
       each final n= n, each final nS=nSL, c=c_l_z[:,1:m,:], V=V_reac, Ndot_l_transfer=Ndot_z[:,1:m,:],  gamma=gamma_z);
  SI.MolarFlowRate Ndot_R[n,m,nSL] = if reaction[1,1].film then reaction.Ndot else  fill( 0,n,m,nSL);
  SI.Energy Qdot_R[n,m] = if reaction[1,1].film then reaction.deltaH_R else fill( 0,n,m);
  SI.Volume V_reac[n,m] "reaction volume";
  Real gamma_z[n,m,nSL]= activityCoeff.gamma;
  SI.Pressure p_film[n,m+1]; //no pressure gradient in film
  MediumLiquid.ActivityCoefficient activityCoeff[n,m](T=T_z[:,1:m],x_l=x_l_z[:,1:m,:]);
  SI.Concentration c_l_z[n,m+1,nSL];

equation
for j in 1:n loop
  for z in 1:m loop
    for i in 1:nSL loop
  diffCoeffLiqMatrix[j,z].eta[i] = eta_comp[j,i];
end for;
end for;
end for;
    for z in 1:m loop
      t_z[:,m+1-z] = t[:]*((z/m)^(1/e)-((z-1)/m)^(1/e));
      t_temp_z[:,m+1-z] = t_temp[:]*((z/m)^(1/e)-((z-1)/m)^(1/e));
    end for;

  /*** link between interface and bulk, vapour film not discretized ***/
  Ndot_v_interface=Ndot_v_transfer;
  Edot_v_interface=Edot_v_transfer;
  Ndot_z[:,1,:] = Ndot_l_interface;
  Ndot_z[:,m+1,:] = Ndot_l_transfer;
  x_l_z[:,1,:] = x_l_star;
  x_l_z[:,m+1,:] = x_l;
  Edot_z[:,1] = Edot_l_interface;
  Edot_z[:,m+1] = Edot_l_transfer;
  T_z[:,1] = T_star;
  T_z[:,m+1] = T_l;
  /*** balance equations ***/
for j in 1:n loop
for z in 1:m loop
   Ndot_z[j,z,:] - Ndot_z[j,z+1,:] + Ndot_R[j,z,:] = zeros(nSL);
   Edot_z[j,z] - Edot_z[j,z+1] + Qdot_R[j,z] = 0;
end for;
end for;

 for j in 1:n loop
    vector_vap[j,:] = (-x_v[j,:] + x_v_star[j,:]);
  end for;

  for z in 1:m loop
    delta_x[:,z,:] = x_l_z[:,z+1,:] - x_l_z[:,z,:];
  end for;

  for j in 1:n loop
    for z in 1:m loop
    B[j,z,1:nSLi,1:nSLi]*Ndot_z[j,z,1:nSLi] =- A_I[j] * (1/propsLiq[j].v) *(thermoFactorL[j,z].Gamma[1:nSLi,1:nSLi]* delta_x[j,z,1:nSLi]/t_z[j,z] + Modelica.Constants.F/(Modelica.Constants.R*T_z[j,z])*c_l_z[j,z,1:nSLi].*MediumLiquid.ic[1:nSLi]*pot_diff_z[j,z]) +sum(Ndot_z[j,z,:])*B[j,z,1:nSLi,1:nSLi]*x_l_z[j,z+1,1:nSLi];
    if sum(abs(MediumLiquid.ic))==0 then
         //no charged particles are in the solution
         pot_diff_z[j,z] =0;
    else
      sum(MediumLiquid.ic[:].*x_l_z[j,z,:])=0; //electroneutrality condition
      pot_diff_z[j,z]= -sum(MediumLiquid.ic[:].*(-c_l_z[j,z+1,:]+c_l_z[j,z,:]))*(Modelica.Constants.R*T_z[j,z])/sum(MediumLiquid.ic[:].*MediumLiquid.ic[:].*c_l[j,:])/Modelica.Constants.F;
    end if;
    sum(x_l_z[j,z,:])=1;
    end for;
  end for;

    // Vapour side
 for j in 1:n loop
     R_v[j,1:nSV-1,1:nSV-1]* Ndot_v_transfer[j,1:nSV-1] =Ndot_v_tot[j] *R_v[j,1:nSV-1,1:nSV-1]*   x_v[j,1:nSV-1] +  A_I[j] * propsVap[j].rho/propsVap[j].MM*thermoFactorV[j].Gamma[1:nSV-1,1:nSV-1]*vector_vap[j,1:nSV-1];
     sum(x_v_star[j,:])=1;
     Ndot_v_tot[j] = sum(Ndot_v_transfer[j,:]);
  end for;

  /*** phase equilibrium ***/
  for j in 1:n loop
   for i in 1:nS loop
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

   /*** energy transfer ***/
   /*** vapour film ***/
   for j in 1:n loop
     for i in 1:nSV loop
    Ndot_fromV[j,i] = -1*min(0,Ndot_v_transfer[j,i]);
    Ndot_toV[j,i] = max(0,Ndot_v_transfer[j,i]);
    x_transfer_fromV[j,i] = Ndot_fromV[j,i]/max(1e-9,sum(Ndot_fromV[j,:]));
    x_transfer_toV[j,i] = Ndot_toV[j,i]/max(1e-9,sum(Ndot_toV[j,:]));
    end for;
   Qdot_v_transfer[j] = alphaVap[j]*A_I[j]*(T_star[j] - T_v[j]);
   Edot_v_transfer[j] = Qdot_v_transfer[j]- sum(Ndot_fromV[j,:])*h_transfer_fromV[j] + sum(Ndot_toV[j,:])*h_transfer_toV[j];
   end for;

   /*** liquid film ***/
   /*** To Do: Energiefluss auf Grund von Molenstrom ***/
   for j in 1:n loop
     for z in 1:m loop
       for i in 1:nSL loop
         Ndot_rightToLeft[j,z,i] = -1*min(0,Ndot_z[j,z,i]);
         Ndot_leftToRight[j,z,i] = max(0,Ndot_z[j,z,i]);
         x_transfer_rightToLeft[j,z,i] = Ndot_rightToLeft[j,z,i]/max(1e-9,sum(Ndot_rightToLeft[j,z,:]));
         x_transfer_leftToRight[j,z,i] = Ndot_leftToRight[j,z,i]/max(1e-9,sum(Ndot_leftToRight[j,z,:]));
    end for;
       Edot_z[j,z] = lambda[j] * A_I[j] * (T_z[j,z]- T_z[j,z+1])/t_temp_z[j,z] + sum(Ndot_z[j,z,:])*h_l_z[j,z] - sum(Ndot_z[j,z+1,:])*h_l_z[j,z+1];
     end for;
     end for;

     /*** reaction ***/
     for j in 1:n loop
     V_reac[j,:] = A_I[j] * t_z[j,:];
     for z in 1:m+1 loop
       for i in 1:nSL loop
     c_l_z[j,z,i] = x_l_z[j,z,i]*mediumLiquid[j,z].v;
       end for;
     end for;
     p_film[j,:] = fill(p_v[j],m+1);
     end for;

 assert(not homotopyMethod.bool_Ndot_inter,"this film model does not support the homotopy method Ndot");
 assert(not homotopyMethod.bool_Edot_inter,"this film model does not support the homotopy method Edot");
end FilmDiscretization;
