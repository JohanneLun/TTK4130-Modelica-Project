within WasteWater.ASM1;
model OCI "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Modelica.Blocks.Interfaces.RealInput AE annotation (Placement(transformation(extent={{-110,84},
            {-90,104}})));
  Modelica.Blocks.Interfaces.RealInput PE annotation (Placement(transformation(extent={{-110,38},
            {-90,58}})));
  Modelica.Blocks.Interfaces.RealInput SP annotation (Placement(transformation(extent={{-110,
            -82},{-90,-102}})));
  Modelica.Blocks.Interfaces.RealInput EC annotation (Placement(transformation(extent={{-110,
            -36},{-90,-56}})));
  Modelica.Blocks.Interfaces.RealInput ME annotation (Placement(transformation(extent={{-110,10},
            {-90,-10}})));
  Modelica.Blocks.Interfaces.RealOutput OCI annotation (Placement(transformation(extent={{92,-10},
            {112,10}})));
equation
  OCI = AE + PE + 5*SP +3*EC + ME;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end OCI;
