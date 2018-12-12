within ThermalSeparation.Reaction.Homogeneous;
model NoKinetic
  "no reaction kinetic, but reaction enthalpy - mass balance optional"
  extends ThermalSeparation.Reaction.BaseReaction;
    parameter Integer nR=1 "number of reactions";
  Units.MolarEnthalpy h_R[nR] "molar reaction enthalpy, negative if exotherm";
  parameter Real a_R[nR] = {38.814} "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R[nR] = {-89155} "h_R [J/mol] = a_R * T [K] + b_R";

  /*** if reaction is included in molar balance ***/
  parameter Boolean moleBalance = false
    "true, if reaction is included in mole balance ";
  parameter Integer reacComp[nR]={1}
    "index of component in the component vector to which the factor refers" annotation(Dialog(enable=moleBalance));
  parameter Real factor[nR,nS] = {{0, -1, -1}}
    "Ndot_reaction[j,i] = factor[i]*Ndot_l_transfer[j,reacComp]" annotation(Dialog(enable=moleBalance));

    SI.Energy deltaH_singleReac[nR];

protected
      SI.MolarFlowRate Ndot_singleReac[
                                 nR,nS];
equation

 if moleBalance then
    for i in 1:nS loop
      for m in 1:nR loop
    Ndot_singleReac[m,i] = factor[m,i]*Ndot_l_transfer[reacComp[nR]];
      end for;
        Ndot[i] = sum(Ndot_singleReac[:,i]);
      end for;
else
    Ndot_singleReac = zeros(nR,nS);
    Ndot = zeros(nS);
end if;

    for m in 1:nR loop
      h_R[m] =a_R[m]*T + b_R[m];
      deltaH_singleReac[m] = - Ndot_l_transfer[reacComp[m]]*h_R[m];
    end for;
    deltaH_R = sum(deltaH_singleReac[:]);

end NoKinetic;
