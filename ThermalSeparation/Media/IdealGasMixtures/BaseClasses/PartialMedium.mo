within ThermalSeparation.Media.IdealGasMixtures.BaseClasses;
partial package PartialMedium
  "Partial medium properties (base package of all media packages)"

  // explicit derivative functions for finite element models

  package Choices "Types, constants to define menu choices"
    package Init
      "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer NoInit=1;
      constant Integer InitialStates=2;
      constant Integer SteadyState=3;
      constant Integer SteadyMass=4;
      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                NoInit "NoInit (no initialization)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                InitialStates "InitialStates (initialize medium states)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                SteadyState "SteadyState (initialize in steady state)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                SteadyMass
              "SteadyMass (initialize density or pressure in steady state)"));
      end Temp;
      annotation (preferedView="text");
    end Init;

    package ReferenceEnthalpy
      "Type, constants and menu choices to define reference enthalpy, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer ZeroAt0K=1;
      constant Integer ZeroAt25C=2;
      constant Integer UserDefined=3;
      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.
                ZeroAt0K
              "The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.
                ZeroAt25C
              "The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.
                UserDefined
              "The user-defined reference enthalpy is used at 293.15 K (25 degC)"));

      end Temp;
      annotation (preferedView="text");
    end ReferenceEnthalpy;

    package ReferenceEntropy
      "Type, constants and menu choices to define reference entropy, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer ZeroAt0K=1;
      constant Integer ZeroAt0C=2;
      constant Integer UserDefined=3;
      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                ZeroAt0K "The entropy is 0 at 0 K (default)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                ZeroAt0C "The entropy is 0 at 0 degC",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                UserDefined
              "The user-defined reference entropy is used at 293.15 K (25 degC)"));

      end Temp;
      annotation (preferedView="text");
    end ReferenceEntropy;

    package pd
      "Type, constants and menu choices to define whether p or d are known, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer default=1;
      constant Integer p_known=2;
      constant Integer d_known=3;

      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.pd.default
              "default (no boundary condition for p or d)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.pd.p_known
              "p_known (pressure p is known)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.pd.d_known
              "d_known (density d is known)"));
      end Temp;
      annotation (preferedView="text");
    end pd;

    package Th
      "Type, constants and menu choices to define whether T or h are known, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer default=1;
      constant Integer T_known=2;
      constant Integer h_known=3;

      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Th.default
              "default (no boundary condition for T or h)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Th.T_known
              "T_known (temperature T is known)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Th.h_known
              "h_known (specific enthalpy h is known)"));
      end Temp;
      annotation (preferedView="text");
    end Th;

    package Explicit
      "Type, constants and menu choices to define the explicitly given state variable inputs"

      constant Integer dT_explicit=0 "explicit in density and temperature";
      constant Integer ph_explicit=1
        "explicit in pressure and specific enthalpy";
      constant Integer ps_explicit=2
        "explicit in pressure and specific entropy";
      constant Integer pT_explicit=3 "explicit in pressure and temperature";

      type Temp
        "Temporary type with choices for menus (until enumerations are available)"
        extends Integer(min=0,max=3);
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.dT_explicit
              "explicit in d and T",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.ph_explicit
              "explicit in p and h",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.ps_explicit
              "explicit in p and s",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.pT_explicit
              "explicit in p and s"));
      end Temp;
    end Explicit;
  end Choices;
  annotation (Documentation(info="<html>
<p>
<b>PartialMedium</b> is a package and contains all <b>declarations</b> for
a medium. This means that constants, models, and functions
are defined that every medium is supposed to support
(some of them are optional). A medium package
inherits from <b>PartialMedium</b> and provides the
equations for the medium. The details of this package
are described in  
<a href=\"Modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.
</p>
</html>
",
 revisions="<html>
  
</html>"));

end PartialMedium;
