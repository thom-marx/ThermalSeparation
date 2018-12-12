within ThermalSeparation.Components.SourcesSinks.Enumerations;
type FlowOption = enumeration(
    FlowVdot "Vdot: variable flow corresponds to volume flow rate in m3/s",
    FlowNdot "Ndot: variable flow corresponds to molar flow rate in mol/s",
    FlowMdot "Mdot: variable flow corresponds to mass flow rate in kg/s")
  "Enumeration with choices for different flow quantities";
