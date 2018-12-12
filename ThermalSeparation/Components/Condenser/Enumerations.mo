within ThermalSeparation.Components.Condenser;
package Enumerations
  type OutletTempOption = enumeration(
      T_in "outlet temperature equals inlet temperature",
      T_set "outlet temperature equals T_set",
      T_subcool "subcooled liquid")
    "Different options for liquid outlet temperature";

end Enumerations;
