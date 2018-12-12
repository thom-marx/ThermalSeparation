within ThermalSeparation.Components.Reboiler;
type InitOption = enumeration(
    initEQ "init with Ndot_transfer = 0",
    initX "init all mole fraction in vapour and liquid",
    initX_Tsat "init with T=T_sat(p)")
  "Enumeration with choices for model structure in distributed pipe model";
