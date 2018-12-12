within ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour;
model KineticTheory
  "theoretical equation for the mutual diffusion coefficient in a low-pressure binary gas mixture (ideal gas!)"
  //z.B. Taylor, p. 68
  extends
    ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour.BaseDiffusionCoeffGas;

protected
  parameter Integer counter1[nS-1]={nS-i+1 for i in 2:nS};//für nS=4: {0,3,2,1};
  Integer counter[nS];
  final parameter Real C= 1.883e-2;
  Real sigma_mix[a] "Lennard-Jones force constant), constant value";
  Real epsilon_k_mix[a]
    "Lennard-Jones force constant (epsilon/k), constant value";
  Real omega[n,a] "collision integral";
  Modelica.Blocks.Tables.CombiTable1D table_omega[n,a](each table=collisionIntegral);

equation
  for j in 1:n loop
   for i in 1:a loop
    table_omega[j,i].u = {T[j]/epsilon_k_mix[i]};
    omega[j,i] = table_omega[j,i].y[1];
   end for;
   end for;

  counter[1]=0;
  counter[2:nS] = counter1;
for j in 1:n loop
for i in 1:nS-1 loop
  for k in i+1:nS loop
    //the factor 1000 is due to unit conversion
 D[j,k-i+ sum(counter[1:i])] = ((MediumVapour.MMX[i] + MediumVapour.MMX[k])/(MediumVapour.MMX[i] * MediumVapour.MMX[k] *1000))^0.5 * C * T[j]^(3/2)/(p[j]*sigma_mix[k-i+ sum(counter[1:i])]^2 * omega[j,k-i+ sum(counter[1:i])]);
  end for;
  end for;
  end for;

 for i in 1:nS-1 loop
  for k in i+1:nS loop
        sigma_mix[k-i+ sum(counter[1:i])] = 0.5*(MediumVapour.sigma[i] + MediumVapour.sigma[k]);
    epsilon_k_mix[k-i+ sum(counter[1:i])]=sqrt(MediumVapour.epsilon_k[i]*MediumVapour.epsilon_k[k]);
      end for;
  end for;

  annotation (Documentation(info="<html>
<p>The diffusion coefficient is calculated by a correlation which is based on the kinetic theory of gases. The equation can be found for instance in [1] and [2]. Values for sigma and epsilon/k for the binary pairs are calculated using the values for sigma and epsilon for the pure substances. These values must be supplied in the medium model. The values for many different substances can be found in [3]. If the values for the pure substances are not available, they can be estimated using the critical data (see [2]).</p>
<p><br/>References:</p>
<p>[1] Taylor, R. and Ross, R.: Multicomponent Mass Transfer, Wiley, 1993, p. 68 </p>
<p>[2] Sherwood, T. et al.: Mass Transfer, McGraw-Hill, 1952, p. 19-23</p>
<p>[3] Svehla, R. A.: Estimated Viscosities and ..., NASA Tech. Rep. R-132, 1962</p>
</html>"));
end KineticTheory;
