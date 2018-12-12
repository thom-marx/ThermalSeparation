within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.TrayColumn;
model Stichlmair
    extends BaseTray;
  //Stichlmair: Grundlagen der Dimensionierung des Gas/Fluesssigkeit-Kontaktapparates Bodenkolonne, S. 109
  parameter Boolean a_init_const = true
    "true if a constant value is assumed for the first second - improves numerical robustness";
  parameter Units.VolumetricArea a_init = 400
    "value for the volumetric area during the first second";
protected
  Units.VolumetricArea a_T[n]
    "specific heat and mass transfer area, Tropfenregime";
    Units.VolumetricArea a_B[n]
    "specific heat and mass transfer area, Blasenregime";
    Units.VolumetricArea a_B_star[n];
    Units.VolumetricArea a_T_star[n];
equation
  for j in 1:n loop
    a_T[j] = F[j]^2/(2*sigma[j]*geometry.phi^2) * (1-(F[j]/F_max[j])^0.28);
    a_T_star[j] = (0.7*F_max[j])^2/(2*sigma[j]*geometry.phi^2)*(1-0.7^0.28);
    a_B[j] = 6*sqrt((rho[j]-rho_v[j])*Modelica.Constants.g_n/(6*sigma[j])) *(F[j]/F_max[j])^0.28;
    a_B_star[j] = 6*sqrt((rho[j]-rho_v[j])*Modelica.Constants.g_n/(6*sigma[j])) * 0.7^0.28;
    if a_init_const then
    if time < 1 then
      a[j]=a_init;
      else
    a[j]= 400;// if F[j]/F_max[j]>0.7 then a_T[j] else a_B[j] - (F[j]/F_max[j]/0.7)^2*(a_B_star[j] - a_T_star[j]);
    end if;
    else
      a[j]= 400;// if F[j]/F_max[j]>0.7 then a_T[j] else a_B[j] - (F[j]/F_max[j]/0.7)^2*(a_B_star[j] - a_T_star[j]);
      end if;
        end for;

//   extends BaseTray;
//   //Stichlmair: Grundlagen der Dimensionierung des Gas/Fluesssigkeit-Kontaktapparates Bodenkolonne, S. 109
//   parameter Boolean a_init_const = true
//     "true if a constant value is assumed for the first second - improves numerical robustness";
//   parameter Units.VolumetricArea a_init = 400
//     "value for the volumetric area during the first second";
// protected
//   Units.VolumetricArea a_T[n]
//     "specific heat and mass transfer area, Tropfenregime";
//     Units.VolumetricArea a_B[n]
//     "specific heat and mass transfer area, Blasenregime";
//     Units.VolumetricArea a_B_star[n];
//     Units.VolumetricArea a_T_star[n];
// equation
//   for j in 1:n loop
//     a_T[j] = F[j]^2/(2*sigma[j]*geometry.phi^2) * (1-(F[j]/F_max[j])^0.28);
//     a_T_star[j] = (0.7*F_max[j])^2/(2*sigma[j]*geometry.phi^2)*(1-0.7^0.28);
//     a_B[j] = 6*sqrt((rho[j]-rho_v[j])*Modelica.Constants.g_n/(6*sigma[j])) *(F[j]/F_max[j])^0.28;
//     a_B_star[j] = 6*sqrt((rho[j]-rho_v[j])*Modelica.Constants.g_n/(6*sigma[j])) * 0.7^0.28;
//     if a_init_const then
//     if time < 1 then
//       a[j]=a_init;
//       else
//     a[j]= if F[j]/F_max[j]>0.7 then a_T[j] else a_B[j] - (F[j]/F_max[j]/0.7)^2*(a_B_star[j] - a_T_star[j]);
//     end if;
//     else
//       a[j]=  if F[j]/F_max[j]>0.7 then a_T[j] else a_B[j] - (F[j]/F_max[j]/0.7)^2*(a_B_star[j] - a_T_star[j]);
//       end if;
//         end for;
end Stichlmair;
