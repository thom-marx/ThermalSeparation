within ThermalSeparation.Components.SourcesSinks.Enumerations;
type SplitOption = enumeration(
    out1_over_in "split = Vdot_out1 / V_in",
    out2_over_in "split = Vdot_out2 / V_in",
    out2_over_out1 "split = Vdot_out2 / Vdot_out1",
    Nout1_over_Nin "split = Ndot_out1 / Ndot_in") "different splitting options";
