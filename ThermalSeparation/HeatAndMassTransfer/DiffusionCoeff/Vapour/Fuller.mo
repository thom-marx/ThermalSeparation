within ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour;
model Fuller "correlation of Fuller"
  extends
    ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour.BaseDiffusionCoeffGas;

  //Werte für N2, H2O, CO2, O2, Taylor S. 69
  //parameter Real V[:]={18.5, 13.1, 16.3, 26.7};
protected
  parameter Integer counter1[nS-1]={nS-i+1 for i in 2:nS};//für nS=4: {0,3,2,1};
  final parameter Real C= 1.013e-2;
  Integer counter[nS];

equation
  counter[1]=0;
  counter[2:nS] = counter1;
for j in 1:n loop
for i in 1:nS-1 loop
  for k in i+1:nS loop
 D[j,k-i+ sum(counter[1:i])] = ((MediumVapour.MMX[i] + MediumVapour.MMX[k])/(MediumVapour.MMX[i] * MediumVapour.MMX[k]*1000))^0.5 * C * T[j]^1.75/p[j]/(MediumVapour.V[i]^(1/3) + MediumVapour.V[k]^(1/3))^2;
  end for;
  end for;
  end for;

  annotation (Documentation(info="<html>
<p>The diffusion coefficient is calculated by a correlation from Fuller et al. [1], [2]. It calcuates the binary diffusion coeff. from the molar masses of the substances, the system pressure and temperature and the molecular diffusion volumes for the substances. The molecular diffusion volumes are calculated by summing the atomic contributions provided in [2].</p>
<p><br/>References:</p>
<p>[1] Fuller et al.: Diffusion of halogenated hydrocarbons in helium. The effect of structure on collision cross sections, J. Phys. Chem., Vol. 73 (1969), pp. 3679-3685</p>
<p>[2] Fuller et al.: A new method for rediction of binary gas-phase diffusion coefficients, Ind. Eng. Chem., Vol 58 (1966), pp. 19-27</p>
</html>"));
end Fuller;
