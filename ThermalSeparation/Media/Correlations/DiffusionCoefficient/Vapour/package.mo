within ThermalSeparation.Media.Correlations.DiffusionCoefficient;
package Vapour

    constant Real collisionIntegral[:,:] = [0.3, 2.662;
                                            0.4, 2.318;
                                            0.5, 2.066;
                                            0.6, 1.877;
                                            0.7, 1.729;
                                             0.8, 1.612;
                                             0.9, 1.517;
                                             1.0, 1.439;
                                             1.1, 1.375;
                                             1.2, 1.320;
                                             1.3, 1.273;
                                             2.0, 1.075;
                                             4.0, 0.8836;
                                             5.0, 0.8422;
                                       400, 0.4170];


  annotation (Documentation(info="<html>
<p>For gases the diffusion coeffient D_12 equals the diffusion coefficient D_21.</p>
</html>"));
end Vapour;
