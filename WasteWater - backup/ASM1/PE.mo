within WasteWater.ASM1;
model PE "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Modelica.Blocks.Interfaces.RealInput Qa annotation (Placement(transformation(extent={{-110,30},
            {-90,50}})));
  Modelica.Blocks.Interfaces.RealInput Qr annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Modelica.Blocks.Interfaces.RealInput Qw annotation (Placement(transformation(extent={{-110,
            -30},{-90,-50}})));
  Modelica.Blocks.Interfaces.RealOutput PE(start=0) annotation (Placement(transformation(extent={{92,-10},
            {112,10}})));
  Real T(start=1e-3);
equation
  der(T) = 1.0;
  der(PE*T) = 0.004 * Qa + 0.008 * Qr + 0.05 * Qw;

end PE;
