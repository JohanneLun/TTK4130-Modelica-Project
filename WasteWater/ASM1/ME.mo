within WasteWater.ASM1;
model ME "Ideal sensor to measure total suspended solids concentration (ASM1)"

  extends WasteWater.Icons.sensor_TSS;
  Modelica.Blocks.Interfaces.RealInput Kla3 annotation (Placement(transformation(extent={{-110,00},
            {-90,-20}})));
  Modelica.Blocks.Interfaces.RealInput Kla4 annotation (Placement(transformation(extent={{-110,-40},
            {-90,-60}})));
  Modelica.Blocks.Interfaces.RealInput Kla5 annotation (Placement(transformation(extent={{-110,-80},
            {-90,-100}})));

  Modelica.Blocks.Interfaces.RealOutput ME(start=0) annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));

  Real From_ax_1;
  Real From_ax_2;
  Real From_ax_3;
  Real T(start=1e-3);
  Real V = 1333;

equation
  der(T) = 1;
  From_ax_1 = if
                (Kla3 < 20) then 0.005*V else 0;
  From_ax_2 = if
                (Kla4 < 20) then 0.005*V else 0;
  From_ax_3 = if (Kla5 < 20) then 0.005*V else 0;

  der(ME*T) = 24*(From_ax_1 + From_ax_2 + From_ax_3);

end ME;
