within ThermalSeparation.Components.Columns.BaseClasses.Initialization;
model Init_FluidProp "init for FluidProp"
extends BaseInit;
initial equation

            for j in 1:n loop
             der(rho_l[j])=0;
             der(rho_v[j])=0;
            for i in 1:nSL loop
             if inertLiquid[i] then
                x_l[j,i] = x_l_start[j,i];
             end if;

          end for;

               for i in 3:nSV loop
              x_v[j,i] = x_v_start[j,i];
          end for;
       sum(x_l[j,:])=1;
       sum(x_v[j,:])=1;
    end for;
     T_v = T_v_start;
     T_l = T_l_start;
    p_v[1:n] = p_v_start;
    for i in 1:nSV loop
      if not inertVapour[i] then
        Ndot_v_transfer[:,i]=zeros(n);
      end if;
      end for;

end Init_FluidProp;
