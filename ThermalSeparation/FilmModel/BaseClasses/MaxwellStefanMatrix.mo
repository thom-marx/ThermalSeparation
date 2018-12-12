within ThermalSeparation.FilmModel.BaseClasses;
model MaxwellStefanMatrix "R or B matrix for Maxwell-Stefan mass transport"
  parameter Integer nS = 2
                       annotation(Dialog(enable=false));

  input SI.MoleFraction x[nS];
  input Real dummy[nS,nS] "dummy variable to be replaced either by k or by D";
  output Real matrix[nS-1,nS-1];

   // If the medium-mixture is considered to be diluted, and the solvent is the last element in the medium vector
   // the R-matrix reduces to a diagonal matrix, where all off-diagonal entries are zero.
   parameter Boolean diluted = false
    "true, if the solution is diluted AND the solvent is the last element in the medium vector";

equation
    for i in 1:nS-1 loop
      for m in 1:nS-1 loop
        if i ==m then
          matrix[i,m] = if diluted then 1/dummy[i,nS] else x[i]/dummy[i,nS] + sum(x[:]./dummy[i,:]) - x[i]/dummy[i,i];
        else
          matrix[i,m] = if diluted then 0 else -x[i]* ( 1/dummy[i,m] - 1/dummy[i,nS]);
        end if;
      end for;
        end for;

end MaxwellStefanMatrix;
