within ThermalSeparation.PhaseEquilibrium.HenrysLaw;
model Exponential "exponential approach for temperature dependency"
  extends BaseHenry;
equation
  for j in 1:n loop
 for i in 1:nSV loop
        He[j,i] = if henry_temp then MediumVapour.henry_H[i] * exp(MediumVapour.henry_C[i]*(1/T[j] - 1/MediumVapour.henry_T_norm)) else MediumVapour.henry_H[i];
 end for;
 end for;
  annotation (Documentation(info="<html>
<p>The Henry coefficient Hi is calculated as a function of the temperature if the parameter henry_temp is true. In this case additionally to a Henry coefficent H_consti for a reference temperature T_norm also the constant Ci is used. Those three values shall be supplied in the vapour medium model.</p>
</html>"));
end Exponential;
