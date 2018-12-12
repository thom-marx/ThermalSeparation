within ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor.ThieleModulus;
model FirstOrder "thiele modulus for first order reaction"
extends BaseThiele;
 SI.DiffusionCoefficient D_eff_cat[nS]
    "effective diffusion coefficient inside the catalyst";
    parameter Real tau "catalyst tortuosity";
    parameter Real epsilon "catalyst porosity";
MediumLiquid.DiffusionCoefficient DiffCoeff(T=propsLiq.T, p=propsLiq.p, eta=propsLiq.eta_comp, x=propsLiq.x);
SI.DiffusionCoefficient D_eff[nS]= DiffCoeff.D_diluted
    "effective diffusion coefficients";
equation
  for i in 1:nS loop
    D_eff_cat[i] = epsilon/tau * D_eff[i];
   thiele[i]=if MediumLiquid.nu[1,i]==0 then 1e6 else V_cat/S_x*sqrt(k_v/D_eff_cat[i]);
  end for;

end FirstOrder;
