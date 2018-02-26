within WasteWater.ASM1;
model Aeration_energy
  "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Modelica.Blocks.Interfaces.RealInput Kla3 annotation (Placement(transformation(extent={{-110,30},
            {-90,50}})));
  Modelica.Blocks.Interfaces.RealInput Kla4 annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Modelica.Blocks.Interfaces.RealInput Kla5 annotation (Placement(transformation(extent={{-110,
            -30},{-90,-50}})));
  Modelica.Blocks.Interfaces.RealOutput AE(start=0) annotation (Placement(transformation(extent={{92,-10},
            {112,10}})));
  Real T(start=1e-3);
equation
  der(T) = 1.0;
  der(AE*T) = 2/1.8/1000*1333*(Kla3 + Kla4 + Kla5);

  connect(AE, AE) annotation (Line(
      points={{102,0},{102,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Documentation(info="This component measures the flow of an ASM1 wastewater stream and provides
the result as output signal (to be further processed with blocks of
the Modelica.Blocks library).
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics));
end Aeration_energy;
