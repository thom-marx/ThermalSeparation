within ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy;
model NoHomotopy
extends ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.BaseHomotopy(
      useHomotopy=false);

Boolean whichHomotopy[6]=fill(false,6);
Boolean bool_Ndot_inter=if not whichHomotopy[1] then false else true
    "true if homotopy is applied on molar flow rate across phase boundary";
Boolean bool_Edot_inter=if not whichHomotopy[2] then false else true
    "true if homotopy is applied on heat flow rate across phase boundary";
Boolean bool_dp=if not whichHomotopy[3] then false else true
    "true if homotopy is applied on pressure loss";
Boolean bool_K=if not whichHomotopy[4] then false else true
    "true is homotopy is applied on equilibrium constant";

Boolean bool_rho=if not whichHomotopy[5] then false else true "true is homotopy is applied on liquid and vapour density";
Boolean bool_h=if not whichHomotopy[6] then false else true "true is homotopy is applied on liquid and vapour enthalpy";

end NoHomotopy;
