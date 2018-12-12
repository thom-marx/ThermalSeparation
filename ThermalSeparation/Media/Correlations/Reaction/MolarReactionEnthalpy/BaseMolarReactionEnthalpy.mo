within ThermalSeparation.Media.Correlations.Reaction.MolarReactionEnthalpy;
partial model BaseMolarReactionEnthalpy
  "base class for molar reaction enthalpy"
  parameter Integer nR(min=1) "number of reactions" annotation(Dialog(enable=false));
  input SI.Temperature T;

  output Units.MolarEnthalpy h_R[nR]
    "molar reaction enthalpy, negative if exotherm";

end BaseMolarReactionEnthalpy;
