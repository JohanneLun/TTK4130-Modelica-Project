within WasteWater.ASM1;
model sensor_Q_INF
  "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  extends Interfaces.stoichiometry;
  Real S_Nkj;
  Real SS_e;
  Real BOD;
  Real COD;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput IQ(start=0) annotation (Placement(transformation(
          extent={{96,48},{116,68}})));
  Real T(start=1e-3);
equation

  In.Q + Out.Q = 0;
  // eventually abs(In.Q) to be shure to have pos. signal

  S_Nkj = In.Snh + In.Snd + In.Xnd + i_xb*(In.Xbh + In.Xba) + i_xp*(In.Xp + In.Xi);
  SS_e = 0.75*(In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp);
  BOD = 0.25*(In.Ss + In.Xs + (1-f_p)*(In.Xbh + In.Xba));
  COD = In.Ss + In.Si + In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp;
  der(T) = 1.0;
  der(IQ*T) =1/1000*(2 * SS_e + 1 * COD + 30 * S_Nkj + 10 * In.Sno + 2 * BOD)  * In.Q;
  //EQ=1/(T*1000)*der(EQ);

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

    annotation (Placement(transformation(
        origin={0,-98},
        extent={{-10,-10},{10,10}},
        rotation=270)),
    Documentation(info="This component measures the flow of an ASM1 wastewater stream and provides
the result as output signal (to be further processed with blocks of
the Modelica.Blocks library).
"));
end sensor_Q_INF;
