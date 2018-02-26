within WasteWater.ASM1;
model EQ "Ideal sensor to measure dissolved oxygen concentration"

  extends Interfaces.stoichiometry;
  extends WasteWater.Icons.sensor_O2;

  Real S_Nkj;
  Real SS_e;
  Real BOD;
  Real COD;

    /*Real Bss;
  Real Bcod;
  Real Bnk;
  Real Bno;
  Real Bbod;*/

  /*Real Q;
  Real Si;
  Real Ss;
  Real Xi;
  Real Xs;
  Real Xbh;
  Real Xba;
  Real Xp;
  Real So;
  Real Sno;
  Real Snh;
  Real Snd;
  Real Xnd;
  Real Salk;*/

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput EQ(start=0) annotation (Placement(transformation(
          extent={{96,48},{116,68}})));
  Real T(start=1e-3);

equation
  /*
  Q = In.Q;
  Si = In.Si;
  Ss = In.Ss;
  Xi = In.Xi;
  Xs = In.Xs;
  Xbh = In.Xbh;
  Xba = In.Xba;
  Xp= In.Xp;
  So = In.So;
  Sno = In.Sno;
  Snh = In.Snh;
  Snd = In.Snd;
  Xnd = In.Xnd;
  Salk = In.Salk;*/

  In.Q + Out.Q = 0;
  In.Salk=0;
  In.So=0;

  //X_X,A???

  S_Nkj = In.Snh + In.Snd + In.Xnd + i_xb*(In.Xbh + In.Xba) + i_xp*(In.Xp + In.Xi);
  SS_e = 0.75*(In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp);
  BOD = 0.25*(In.Ss + In.Xs + (1-f_p)*(In.Xbh + In.Xba));
  COD = In.Ss + In.Si + In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp;

  /*Bss = 2;
  Bcod = 1;
  Bnk = 30;
  Bno = 10;
  Bbod = 2;*/

  der(T) = 1.0;
  der(EQ*T) = 1/1000 * (2 * SS_e + 1 * COD + 30 * S_Nkj + 10 * In.Sno + 2 * BOD)  * In.Q;
  //der(EQ*T) = 1/1000 * (Bss * SS_e + Bcod * COD + Bnk * S_Nkj + Bno * In.Sno + Bbod * BOD)  * In.Q;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end EQ;
