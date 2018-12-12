within ThermalSeparation.FilmModel.TrayColumn;
model EquilibriumPI "equilibrium on each stage using a PI controller"
  import ThermalSeparation;
  extends BaseNonEqType(nn=n); //for replaceability
  extends ThermalSeparation.FilmModel.BaseClasses.EquilibriumPI(
                                                              redeclare record
      BaseGeometry =
        ThermalSeparation.Geometry);

  parameter Real eta_murphree[n,nSV] = 0.8*ones(n,nSV) "Murphree efficientcy"          annotation(Dialog(enable=EQ,tab="Heat and Mass Transfer", group = "Equilibrium"));

           Geometry geometry;
equation
    /***for equilibrium model of spray column: deviation from equilibrium calculated using HETP values***/
A_I=fill(400*geometry.H*geometry.A/n,n);
     /***for equilibrium model of tray column: deviation from equilibrium calculated using Murphree efficiency***/
   //first stage
  for i in 1:nSV loop
    x_v_star[1,i] = eta_murphree[1,i] * (x_v_star_eq[1,i] - x_v_in[i]) + x_v_in[i];
  end for;
  //all other stages
  for j in 2:n loop
    for i in 1:nSV loop
       x_v_star[j,i] =  eta_murphree[j,i] * (x_v_star_eq[j,i] - x_v_star_eq[j-1,i]) + x_v_star_eq[j-1,i];
    end for;
    end for;
end EquilibriumPI;
