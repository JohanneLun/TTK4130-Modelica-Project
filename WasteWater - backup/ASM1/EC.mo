within WasteWater.ASM1;
model EC "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Modelica.Blocks.Interfaces.RealInput COD1 annotation (Placement(transformation(extent={{-110,84},
            {-90,104}})));
  Modelica.Blocks.Interfaces.RealInput COD2 annotation (Placement(transformation(extent={{-110,34},
            {-90,54}})));
  Modelica.Blocks.Interfaces.RealInput COD3 annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Modelica.Blocks.Interfaces.RealInput COD4 annotation (Placement(transformation(extent={{-110,
            -54},{-90,-34}})));
  Modelica.Blocks.Interfaces.RealInput COD5 annotation (Placement(transformation(extent={{-110,
            -78},{-90,-98}})));
  Modelica.Blocks.Interfaces.RealOutput EC(start=0) annotation (Placement(transformation(extent={{92,-10},
            {112,10}})));
  Real T(start=1e-3);
equation
  der(T) = 1.0;
  der(EC*T) = 400000/1000 * (COD1 + COD2 + COD3 + COD4 + COD5);

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end EC;
