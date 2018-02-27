within WasteWater.ASM1;
model sensor_SP
  "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-108,
            -52},{-88,-32}})));
  Interfaces.WWFlowAsm1in u annotation (Placement(transformation(extent={{-112,42},
            {-92,62}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput SP(start=0) annotation (Placement(transformation(
        origin={0,-98},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Real T(start=1e-3);
equation
  der(T) = 1;
  In.Q + Out.Q = 0;
  u.Q=0;
  //Q = In.Q;
  // eventually abs(In.Q) to be shure to have pos. signal

  der(SP*T) = 0.75 * (In.Xs + (u.Xi +In.Xi) + In.Xbh + In.Xba)*In.Q;

  In.Si = Out.Si;
  In.Ss = Out.Ss;
  In.Xi = Out.Xi;
  In.Xs = Out.Xs;
  In.Xbh = Out.Xbh;
  In.Xba = Out.Xba;
  In.Xp = Out.Xp;
  In.So = Out.So;
  In.Sno = Out.Sno;
  In.Snh = Out.Snh;
  In.Snd = Out.Snd;
  In.Xnd = Out.Xnd;
  In.Salk = Out.Salk;

  annotation (
    Documentation(info="This component measures the flow of an ASM1 wastewater stream and provides
the result as output signal (to be further processed with blocks of
the Modelica.Blocks library).
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics));
end sensor_SP;
