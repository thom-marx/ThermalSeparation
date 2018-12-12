within ThermalSeparation.BalanceEquations.Base;
partial model BaseBalanceEquations

  replaceable model HomotopyMethod =
      ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy
                       constrainedby
    ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.BaseHomotopy;

end BaseBalanceEquations;
