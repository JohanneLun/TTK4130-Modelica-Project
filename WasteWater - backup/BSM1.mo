within WasteWater;
package BSM1 "Component models for the Activated Sludge Model No.1"
extends Modelica.Icons.Library;

model deni "ASM1 denitrification tank"
  //denitrification tank based on the ASM1 model

  extends WasteWater.Icons.deni;
  extends Interfaces.ASM1base;

  // tank specific parameters
  parameter Modelica.SIunits.Volume V=1000 "Volume of denitrification tank";

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
          extent={{50,40},{60,50}})));
  Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
          extent={{-110,30},{-90,50}})));
  Modelica.Blocks.Interfaces.RealOutput Kla annotation (Placement(
       transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,54})));
  //Modelica.Blocks.Interfaces.RealOutput
equation
  Kla=240;
  aeration = 0;
  // no aeration in this tank //

  // volume dependent dilution term of each concentration

  inputSi = (In.Si - Si)*In.Q/V;
  inputSs = (In.Ss - Ss)*In.Q/V;
  inputXi = (In.Xi - Xi)*In.Q/V;
  inputXs = (In.Xs - Xs)*In.Q/V;
  inputXbh = (In.Xbh - Xbh)*In.Q/V;
  inputXba = (In.Xba - Xba)*In.Q/V;
  inputXp = (In.Xp - Xp)*In.Q/V;
  inputSo = (In.So - So)*In.Q/V;
  inputSno = (In.Sno - Sno)*In.Q/V;
  inputSnh = (In.Snh - Snh)*In.Q/V;
  inputSnd = (In.Snd - Snd)*In.Q/V;
  inputXnd = (In.Xnd - Xnd)*In.Q/V;
  inputSalk = (In.Salk - Salk)*In.Q/V;

  annotation (
    Documentation(info="This component models the ASM1 processes and reactions taking place in an unaerated (denitrification) tank
of a wastewater treatment plant.

The InPort signal gives the tank temperature to the model.

Parameters:

    - all stoichiometric and kinetic parameters of the activated sludge model No.1 (ASM1)
  V - volume of the tank [m3]
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics));
end deni;

model nitri "ASM1 nitrification tank"
  // nitrification (aerated) tank, based on the ASM1 model

  extends WasteWater.Icons.nitri;
  extends Interfaces.ASM1base;

  // tank specific parameters
  parameter Modelica.SIunits.Volume V=1333 "Volume of nitrification tank";

  // aeration system dependent parameters
  parameter Real alpha=0.7 "Oxygen transfer factor";
  parameter Modelica.SIunits.Length de=4.5 "depth of aeration";
  parameter Real R_air=23.5 "specific oxygen feed factor [gO2/(m^3*m)]";
  WWU.MassConcentration So_sat = 8.0 "Dissolved oxygen saturation";

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,-10},
            {110,10}})));
  Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
          extent={{50,40},{60,50}})));
  Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
          extent={{-110,30},{-90,50}})));
  Modelica.Blocks.Interfaces.RealOutput Kla annotation (Placement(
       transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={2,30})));
  Interfaces.AirFlow AirIn annotation (Placement(transformation(extent={{-5,
            -103},{5,-93}})));
equation
  Kla=240;
  //So=2;
  // Temperature dependent oxygen saturation by Simba

  // extends the Oxygen differential equation by an aeration term
  // aeration [mgO2/l]; AirIn.Q_air needs to be in
  // Simulationtimeunit [m3*day^-1]
  //aeration = (alpha*(So_sat - So)/So_sat*de)/V;
  aeration = Kla * (So_sat - So);

  // volume dependent dilution term of each concentration

  inputSi = (In.Si - Si)*In.Q/V;
  inputSs = (In.Ss - Ss)*In.Q/V;
  inputXi = (In.Xi - Xi)*In.Q/V;
  inputXs = (In.Xs - Xs)*In.Q/V;
  inputXbh = (In.Xbh - Xbh)*In.Q/V;
  inputXba = (In.Xba - Xba)*In.Q/V;
  inputXp = (In.Xp - Xp)*In.Q/V;
  inputSo = (In.So - So)*In.Q/V;
  inputSno = (In.Sno - Sno)*In.Q/V;
  inputSnh = (In.Snh - Snh)*In.Q/V;
  inputSnd = (In.Snd - Snd)*In.Q/V;
  inputXnd = (In.Xnd - Xnd)*In.Q/V;
  inputSalk = (In.Salk - Salk)*In.Q/V;

  annotation (
    Documentation(info="This component models the ASM1 processes and reactions taking place in an aerated (nitrification) tank
of a wastewater treatment plant.

The InPort signal gives the tank temperature to the model.

Parameters:

        - all soichiometric and kinetic parameters of the activated sludge model No.1 (ASM1)
  V     - volume of the tank [m3]
  alpha - oxygen transfer factor
  de    - depth of the aeration system [m]
  R_air - specific oxygen feed factor [g O2/(m3*m)]
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics));
end nitri;

  model nitri_tank5 "ASM1 nitrification tank"
    // nitrification (aerated) tank, based on the ASM1 model

    extends WasteWater.Icons.nitri;
    extends Interfaces.ASM1base;

    // tank specific parameters
    parameter Modelica.SIunits.Volume V=1333 "Volume of nitrification tank";

    // aeration system dependent parameters
    parameter Real alpha=0.7 "Oxygen transfer factor";
    parameter Modelica.SIunits.Length de=4.5 "depth of aeration";
    parameter Real R_air=23.5 "specific oxygen feed factor [gO2/(m^3*m)]";
    WWU.MassConcentration So_sat = 8.0 "Dissolved oxygen saturation";

    Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
              -10},{-90,10}})));
    Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
              -10},{110,10}})));
    Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
            extent={{50,40},{60,50}})));
    Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
            extent={{-110,30},{-90,50}})));
    Modelica.Blocks.Interfaces.RealOutput Kla annotation (Placement(
          transformation(extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,30})));

    Interfaces.AirFlow AirIn annotation (Placement(transformation(extent={{-5,
              -103},{5,-93}})));
  equation
    Kla=240;
    //So=2;
    // Temperature dependent oxygen saturation by Simba

    // extends the Oxygen differential equation by an aeration term
    // aeration [mgO2/l]; AirIn.Q_air needs to be in
    // Simulationtimeunit [m3*day^-1]
    //aeration = (alpha*(So_sat - So)/So_sat*de)/V;
    aeration = Kla * (So_sat - So);

    // volume dependent dilution term of each concentration

    inputSi = (In.Si - Si)*In.Q/V;
    inputSs = (In.Ss - Ss)*In.Q/V;
    inputXi = (In.Xi - Xi)*In.Q/V;
    inputXs = (In.Xs - Xs)*In.Q/V;
    inputXbh = (In.Xbh - Xbh)*In.Q/V;
    inputXba = (In.Xba - Xba)*In.Q/V;
    inputXp = (In.Xp - Xp)*In.Q/V;
    inputSo = (In.So - So)*In.Q/V;
    inputSno = (In.Sno - Sno)*In.Q/V;
    inputSnh = (In.Snh - Snh)*In.Q/V;
    inputSnd = (In.Snd - Snd)*In.Q/V;
    inputXnd = (In.Xnd - Xnd)*In.Q/V;
    inputSalk = (In.Salk - Salk)*In.Q/V;

    annotation (
      Documentation(info="This component models the ASM1 processes and reactions taking place in an aerated (nitrification) tank
of a wastewater treatment plant.

The InPort signal gives the tank temperature to the model.

Parameters:

        - all soichiometric and kinetic parameters of the activated sludge model No.1 (ASM1)
  V     - volume of the tank [m3]
  alpha - oxygen transfer factor
  de    - depth of the aeration system [m]
  R_air - specific oxygen feed factor [g O2/(m3*m)]
"),   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
  end nitri_tank5;

model SecClarModTakacs "Secondary Clarifier ASM1 Model based on Takacs"

  extends WasteWater.Icons.SecClar;
  extends SecClar.Takacs.Interfaces.ratios;
  package SCP = SecClar.Takacs;
  import SI = Modelica.SIunits;
  package WI = WasteWater.BSM1.Interfaces;
  package WWU = WasteWater.WasteWaterUnits;
  parameter SI.Length hsc=4.0 "height of secondary clarifier";
  parameter Integer n=10 "number of layers of SC model";
  parameter SI.Length zm=hsc/(1.0*n) "height of m-th secondary clarifier layer";
  parameter SI.Area Asc=1500.0 "area of secondary clarifier";
  parameter WWU.MassConcentration Xt=3000.0 "threshold for X";

  // total sludge concentration in clarifier feed
  WWU.MassConcentration Xf;

  WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},{
            -90,24}})));
  WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,47},
            {112,67}})));
  WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
            -106},{-20,-86}})));
  WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,-106},
            {40,-86}})));

  // layers 1 to 10
  SCP.bottom_layer S1(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    rXs=rXs,
    rXbh=rXbh,
    rXba=rXba,
    rXp=rXp,
    rXi=rXi,
    rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
  SCP.lower_layer S2(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
  SCP.lower_layer S3(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
  SCP.lower_layer S4(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
  SCP.lower_layer S5(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
  SCP.feed_layer S6(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
  SCP.upper_layer S7(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    Xt=Xt) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
  SCP.upper_layer S8(
    zm=zm,
    Asc=Asc,
    Xt=Xt,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
  SCP.upper_layer S9(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    Xt=Xt) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
  SCP.top_layer S10(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    Xt=Xt,
    rXs=rXs,
    rXbh=rXbh,
    rXba=rXba,
    rXp=rXp,
    rXi=rXi,
    rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
equation

  connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
          -2.22045e-15,-74}}));
  connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
          -2.22045e-15,-55}}));
  connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
          -2.22045e-15,-36}}));
  connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
          -2.22045e-15,2}}));
  connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
          -2.22045e-15,21}}));
  connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
          -2.22045e-15,40}}));
  connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
          -2.22045e-15,78}}));
  connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
          -2.22045e-15,-17}}));
  connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
          -2.22045e-15,59}}));
  connect(Feed, S6.In) annotation (Line(points={{-100,14},{-67.5,14},{-67.5,9.8},
          {-35,9.8}}));
  connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-96},{30,-96}}));
  connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
          67.5,57},{102,57}}));
  connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-96},{-30,-96}}));

  // total sludge concentration in clarifier feed
  Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

  // ratios of solid components
  rXs = Feed.Xs/Xf;
  rXbh = Feed.Xbh/Xf;
  rXba = Feed.Xba/Xf;
  rXp = Feed.Xp/Xf;
  rXi = Feed.Xi/Xf;
  rXnd = Feed.Xnd/Xf;

  annotation (
    Documentation(info="This component models an ASM1 10 - layer secondary clarifier model with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Takac`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  Xt  -  threshold value for Xtss [mg/l]

"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics));
end SecClarModTakacs;

model blower "Blower for the aeration of the nitrification tanks"

  extends WasteWater.Icons.blower;
  package WWU = WasteWater.WasteWaterUnits;

  parameter WWU.VolumeFlowRate Q_max=20000 "maximum blower capacity";
  parameter WWU.VolumeFlowRate Q_min=0.0 "minimum blower capacity";
  Real H;
  // this is just a help variable to reduce expressions

  Interfaces.AirFlow AirOut annotation (Placement(transformation(extent={{-20,
            90},{0,110}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(
        origin={98,-30},
        extent={{-10,-10},{10,10}},
        rotation=180)));
equation

  H =0.5*(-Q_min + Q_max) + u*0.5*(-Q_min + Q_max) + Q_min;
  AirOut.Q_air = -(if H > Q_max then Q_max else if H < Q_min then Q_min else H);

  annotation (
    Documentation(info="This component models a blower of a wastewater treatment plant which generates an airflow that is needed
for the nitrification.
The blower is connected to the nitrification tank.
The airflow is controlled by a signal u (-1 <= u <= 1).

Parameter:

  Qmax - maximum blower capacity [m3 Air/d], this is produced when the control signal u is 1 or greater.
  Qmin - minimum blower capacity [m3 Air/d], this is produced when the control signal u is -1 or below.

"));
end blower;

model pump "ASM1 wastewater pump"

  extends WasteWater.Icons.pump;
  package WWU = WasteWater.WasteWaterUnits;

  parameter WWU.VolumeFlowRate Q_min=0.0 "minimum pump capacity";
  parameter WWU.VolumeFlowRate Q_max=20000 "maximum pump capacity";
  Real H;
  // this is just a help variable to reduce expressions

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -43},{-90,-23}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            18},{110,38}})));
  Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
          extent={{-99,15},{-79,35}})));
equation

  //H =0.5*(-Q_min + Q_max) + u*0.5*(-Q_min + Q_max) + Q_min;
  H=u;
  Out.Q = -(if H > Q_max then Q_max else if H < Q_min then Q_min else H);
  //Out.Q = -(if H > Q_max then 14 else if H < Q_min then 11 else H);

  Out.Q + In.Q = 0;
  Out.Si = In.Si;
  Out.Ss = In.Ss;
  Out.Xi = In.Xi;
  Out.Xs = In.Xs;
  Out.Xbh = In.Xbh;
  Out.Xba = In.Xba;
  Out.Xp = In.Xp;
  Out.So = In.So;
  Out.Sno = In.Sno;
  Out.Snh = In.Snh;
  Out.Snd = In.Snd;
  Out.Xnd = In.Xnd;
  Out.Salk = In.Salk;

  annotation (
    Documentation(info="This component models an ASM1 wastewater pump. It generates a wastewater flow
that is controlled by the signal u (-1 <= u <=1).

Parameter:

  Qmax - maximum pump capacity [m3/d], this is produced when the control signal u is 1 or greater.
  Qmin - minimum pump capacity [m3/d], this is produced when the control signal u is -1 or below.

"));
end pump;

model FlowSource "Flowsource"

  extends WasteWater.Icons.FlowSource;
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{88,
            -80},{108,-60}})));
  Modelica.Blocks.Interfaces.RealInput data
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
equation

  Out.Q =-data;

  annotation (
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Ellipse(
          extent={{-54,54},{56,-54}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-8,-54},{-14,-52},{-24,-48},{-32,-44},{-36,-40},{-42,-34},{
              -48,-26},{-50,-20},{52,-20},{50,-26},{46,-32},{42,-36},{38,-40},{
              34,-44},{30,-46},{26,-48},{22,-50},{16,-52},{10,-54},{4,-54},{0,
              -54},{-8,-54}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,95,191},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-54,54},{56,-54}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Rectangle(
          extent={{-4,-52},{4,-74}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,95,191},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{4,-74},{88,-68}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,95,191},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-4,-54},{-4,-74},{88,-74}},
          thickness=0.5),
        Line(
          points={{4,-54},{4,-68},{88,-68}},
          thickness=0.5)}),
    Documentation(info="This component is used to feed an ASM1 wwtp model with flow data from measurement
when e.g. concentration is measured after the primary clarifier.

The dimension of InPort is 1.

  1 volumeflowrate Q of incoming wastewater [m3/d]"));
end FlowSource;

model WWSource "Wastewater source"

  extends WasteWater.Icons.WWSource;
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{88,
            -80},{108,-60}})));
  Modelica.Blocks.Interfaces.RealInput data[14]
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
equation

  Out.Q =-data[1];
  Out.Si =data[2];
  Out.Ss =data[3];
  Out.Xi =data[4];
  Out.Xs =data[5];
  Out.Xbh =data[6];
  Out.Xba =data[7];
  Out.Xp =data[8];
  Out.So =data[9];
  Out.Sno =data[10];
  Out.Snh =data[11];
  Out.Snd =data[12];
  Out.Xnd =data[13];
  Out.Salk =data[14];

  /*Out.Q = -18446;
  Out.Si =30.00;
  Out.Ss =69.50;
  Out.Xi =51.20;
  Out.Xs =202.32;
  Out.Xbh =28.17;
  Out.Xba =0.00;
  Out.Xp = 0.00;
  Out.So = 0.00;
  Out.Sno = 0.00;
  Out.Snh =31.56;
  Out.Snd =6.95;
  Out.Xnd =10.59;
  Out.Salk =7.00;*/

  annotation (
    Documentation(info="This component provides all ASM1 data at the influent of a wastewater treatment plant.
The dimension of InPort is 14.

  1 volumeflowrate Q of incoming wastewater [m3/d]
  2 Si  [g COD/m3]
  3 Ss  [g COD/m3]
  4 Xi  [g COD/m3]
  5 Xs  [g COD/m3]
  6 Xbh [g COD/m3]
  7 Xba [g COD/m3]
  8 Xp  [g COD/m3]
  9 So  [g O2/m3]
 10 Sno [g N/m3]
 11 Snh [g N/m3]
 12 Snd [g N/m3]
 13 Xnd [g N/m3]
 14 Salk[mmol/l]
"));
end WWSource;

model Source_constant "Wastewater source"

  extends WasteWater.Icons.WWSource;
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{88,
            -80},{108,-60}})));
  /*Modelica.Blocks.Interfaces.RealInput data[14]
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));*/
equation

  /*Out.Q =-data[1];
  Out.Si =data[2];
  Out.Ss =data[3];
  Out.Xi =data[4];
  Out.Xs =data[5];
  Out.Xbh =data[6];
  Out.Xba =data[7];
  Out.Xp =data[8];
  Out.So =data[9];
  Out.Sno =data[10];
  Out.Snh =data[11];
  Out.Snd =data[12];
  Out.Xnd =data[13];
  Out.Salk =data[14];*/

  Out.Q = -18446;
  Out.Si =30.00;
  Out.Ss =69.50;
  Out.Xi =51.20;
  Out.Xs =202.32;
  Out.Xbh =28.17;
  Out.Xba =0.00;
  Out.Xp = 0.00;
  Out.So = 0.00;
  Out.Sno = 0.00;
  Out.Snh =31.56;
  Out.Snd =6.95;
  Out.Xnd =10.59;
  Out.Salk =7.00;

  annotation (
    Documentation(info="This component provides all ASM1 data at the influent of a wastewater treatment plant.
The dimension of InPort is 14.

  1 volumeflowrate Q of incoming wastewater [m3/d]
  2 Si  [g COD/m3]
  3 Ss  [g COD/m3]
  4 Xi  [g COD/m3]
  5 Xs  [g COD/m3]
  6 Xbh [g COD/m3]
  7 Xba [g COD/m3]
  8 Xp  [g COD/m3]
  9 So  [g O2/m3]
 10 Sno [g N/m3]
 11 Snh [g N/m3]
 12 Snd [g N/m3]
 13 Xnd [g N/m3]
 14 Salk[mmol/l]
"));
end Source_constant;

model EffluentSink "Receiving water (river)"
  // only for graphical termination in diagram layer, no equation needed

  extends WasteWater.Icons.EffluentSink;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            10},{-90,30}})));
  annotation (
    Documentation(info="This component terminates an ASM1 wastewater treatment plant model e.g. the wastewater flow to the receiving water.
"));
end EffluentSink;

model SludgeSink "Wastesludge sink"
  // only for graphical termination in diagram layer, no equation needed

  extends WasteWater.Icons.SludgeSink;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -22},{-90,-2}})));
  annotation (
    Documentation(info="This component terminates the waste sludge stream of an ASM1 wastewater treatment plant model.
Storage or further sludge treatment is not jet considered."));
end SludgeSink;

model ControlledDivider2 "Controlled flow divider"
  // divides one flow of wastewater into 2 Flows controlled by the

  // input signal u; u=1 means Out1.Q=In.Q and u=0 means Out2.Q=In.Q

  extends WasteWater.Icons.ControlledDivider2;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-111,
            -7},{-91,13}})));
  Interfaces.WWFlowAsm1out Out1 annotation (Placement(transformation(extent={{
            90,15},{110,35}})));
  Interfaces.WWFlowAsm1out Out2 annotation (Placement(transformation(extent={{
            90,-25},{110,-5}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(
        origin={0,-60},
        extent={{-10,-10},{10,10}},
        rotation=90)));
equation

  Out1.Q =-In.Q*u;
  Out2.Q =-In.Q*(1 - u);

  Out1.Si = In.Si;
  Out1.Ss = In.Ss;
  Out1.Xi = In.Xi;
  Out1.Xs = In.Xs;
  Out1.Xbh = In.Xbh;
  Out1.Xba = In.Xba;
  Out1.Xp = In.Xp;
  Out1.So = In.So;
  Out1.Sno = In.Sno;
  Out1.Snh = In.Snh;
  Out1.Snd = In.Snd;
  Out1.Xnd = In.Xnd;
  Out1.Salk = In.Salk;

  Out2.Si = In.Si;
  Out2.Ss = In.Ss;
  Out2.Xi = In.Xi;
  Out2.Xs = In.Xs;
  Out2.Xbh = In.Xbh;
  Out2.Xba = In.Xba;
  Out2.Xp = In.Xp;
  Out2.So = In.So;
  Out2.Sno = In.Sno;
  Out2.Snh = In.Snh;
  Out2.Snd = In.Snd;
  Out2.Xnd = In.Xnd;
  Out2.Salk = In.Salk;

  annotation (
    Documentation(info="This component divides one wastewater flow (ASM1) into two flows which are controlled by the signal u (0...1).
Is u.signal=1, the flow goes to output 1 (Out1) and is u.signal=0, the flow goes to output 2 (Out2).
The concentrations of the outport-flows are equal to the concentration at inport."));
end ControlledDivider2;

model divider2 "Flowdivider"

    // divides one flow of wastewater into 2 Flows; one amount needs to be specified

  extends WasteWater.Icons.divider2;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -7},{-90,13}})));
  Interfaces.WWFlowAsm1out Out1 annotation (Placement(transformation(extent={{
            90,16},{110,36}})));
  Interfaces.WWFlowAsm1out Out2 annotation (Placement(transformation(extent={{
            90,-25},{110,-5}})));
equation

  In.Q + Out1.Q + Out2.Q = 0;

  Out1.Si = In.Si;
  Out1.Ss = In.Ss;
  Out1.Xi = In.Xi;
  Out1.Xs = In.Xs;
  Out1.Xbh = In.Xbh;
  Out1.Xba = In.Xba;
  Out1.Xp = In.Xp;
  Out1.So = In.So;
  Out1.Sno = In.Sno;
  Out1.Snh = In.Snh;
  Out1.Snd = In.Snd;
  Out1.Xnd = In.Xnd;
  Out1.Salk = In.Salk;

  Out2.Si = In.Si;
  Out2.Ss = In.Ss;
  Out2.Xi = In.Xi;
  Out2.Xs = In.Xs;
  Out2.Xbh = In.Xbh;
  Out2.Xba = In.Xba;
  Out2.Xp = In.Xp;
  Out2.So = In.So;
  Out2.Sno = In.Sno;
  Out2.Snh = In.Snh;
  Out2.Snd = In.Snd;
  Out2.Xnd = In.Xnd;
  Out2.Salk = In.Salk;

  annotation (
    Documentation(info=
          "This component divides one ASM1 wastewater flow into two ASM1 wastewater flows."));
end divider2;

model mixer2 "Mixer of two ASM1 characterised flows"

  extends WasteWater.Icons.mixer2;
  Interfaces.WWFlowAsm1in In1 annotation (Placement(transformation(extent={{
            -110,15},{-90,35}})));
  Interfaces.WWFlowAsm1in In2 annotation (Placement(transformation(extent={{
            -110,-25},{-90,-5}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -5},{110,15}})));
equation

  In1.Q + In2.Q + Out.Q = 0;
  Out.Si = (In1.Si*In1.Q + In2.Si*In2.Q)/(In1.Q + In2.Q);
  Out.Ss = (In1.Ss*In1.Q + In2.Ss*In2.Q)/(In1.Q + In2.Q);
  Out.Xi = (In1.Xi*In1.Q + In2.Xi*In2.Q)/(In1.Q + In2.Q);
  Out.Xs = (In1.Xs*In1.Q + In2.Xs*In2.Q)/(In1.Q + In2.Q);
  Out.Xbh = (In1.Xbh*In1.Q + In2.Xbh*In2.Q)/(In1.Q + In2.Q);
  Out.Xba = (In1.Xba*In1.Q + In2.Xba*In2.Q)/(In1.Q + In2.Q);
  Out.Xp = (In1.Xp*In1.Q + In2.Xp*In2.Q)/(In1.Q + In2.Q);
  Out.So = (In1.So*In1.Q + In2.So*In2.Q)/(In1.Q + In2.Q);
  Out.Sno = (In1.Sno*In1.Q + In2.Sno*In2.Q)/(In1.Q + In2.Q);
  Out.Snh = (In1.Snh*In1.Q + In2.Snh*In2.Q)/(In1.Q + In2.Q);
  Out.Snd = (In1.Snd*In1.Q + In2.Snd*In2.Q)/(In1.Q + In2.Q);
  Out.Xnd = (In1.Xnd*In1.Q + In2.Xnd*In2.Q)/(In1.Q + In2.Q);
  Out.Salk = (In1.Salk*In1.Q + In2.Salk*In2.Q)/(In1.Q + In2.Q);

  annotation (
    Documentation(info=
          "This component mixes two flows of wastewater (ASM1) of different concentration and different amount."));
end mixer2;

model mixer3 "Mixer of 3 ASM1 characterised flows"

  extends WasteWater.Icons.mixer3;
  Interfaces.WWFlowAsm1in In1 annotation (Placement(transformation(extent={{
            -110,25},{-90,45}})));
  Interfaces.WWFlowAsm1in In2 annotation (Placement(transformation(extent={{
            -110,-15},{-90,5}})));
  Interfaces.WWFlowAsm1in In3 annotation (Placement(transformation(extent={{
            -110,-55},{-90,-35}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -14},{110,6}})));
equation

  In1.Q + In2.Q + In3.Q + Out.Q = 0;
  Out.Si = (In1.Si*In1.Q + In2.Si*In2.Q + In3.Si*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Ss = (In1.Ss*In1.Q + In2.Ss*In2.Q + In3.Ss*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Xi = (In1.Xi*In1.Q + In2.Xi*In2.Q + In3.Xi*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Xs = (In1.Xs*In1.Q + In2.Xs*In2.Q + In3.Xs*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Xbh = (In1.Xbh*In1.Q + In2.Xbh*In2.Q + In3.Xbh*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Xba = (In1.Xba*In1.Q + In2.Xba*In2.Q + In3.Xba*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Xp = (In1.Xp*In1.Q + In2.Xp*In2.Q + In3.Xp*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.So = (In1.So*In1.Q + In2.So*In2.Q + In3.So*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Sno = (In1.Sno*In1.Q + In2.Sno*In2.Q + In3.Sno*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Snh = (In1.Snh*In1.Q + In2.Snh*In2.Q + In3.Snh*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Snd = (In1.Snd*In1.Q + In2.Snd*In2.Q + In3.Snd*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Xnd = (In1.Xnd*In1.Q + In2.Xnd*In2.Q + In3.Xnd*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Salk = (In1.Salk*In1.Q + In2.Salk*In2.Q + In3.Salk*In3.Q)/(In1.Q + In2.Q
     + In3.Q);

  annotation (
    Documentation(info=
          "This component mixes 3 flows of wastewater (ASM1) of different concentration and different amount."));
end mixer3;

model sensor_COD "Ideal sensor to measure chemical oxygen demand (COD)"

  extends WasteWater.Icons.sensor_COD;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput COD annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0.0;
  COD = In.Si + In.Ss + In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;

  annotation (
    Documentation(info="This component measures the chemical oxygen demand (COD) concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_COD;

model sensor_NH "Ideal sensor to measure ammonium nitrogen"

  extends WasteWater.Icons.sensor_NH;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput Snh annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0;
  Snh = In.Snh;

  annotation (
    Documentation(info="This component measures the ammonium nitrogen concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_NH;

model sensor_NO "Ideal sensor to measure nitrate nitrogen"

  extends WasteWater.Icons.sensor_NO;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput Sno annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0;
  Sno = In.Sno;

  annotation (
    Documentation(info="This component measures the nitrate nitrogen concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_NO;

model sensor_O2 "Ideal sensor to measure dissolved oxygen concentration"

  extends WasteWater.Icons.sensor_O2;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput So annotation (Placement(transformation(
          extent={{88,-10},{108,10}})));
equation

  In.Q = 0;
  So = In.So;

  annotation (
    Documentation(info="This component measures the dissolved oxygen concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"), Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Ellipse(
          extent={{-50,50},{50,-50}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={223,223,159},
          fillPattern=FillPattern.Solid),
        Line(
          points={{0,50},{0,38}},
          thickness=0.5),
        Line(
          points={{-50,0},{38,0}},
          thickness=0.5),
        Line(
          points={{50,0},{38,0}},
          thickness=0.5),
        Line(
          points={{-36,34},{-28,26}},
          thickness=0.5),
        Line(
          points={{34,36},{26,28}},
          thickness=0.5),
        Line(
          points={{0,0},{26,28}},
          thickness=0.5),
        Polygon(
          points={{30,32},{10,24},{24,12},{30,32}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(extent={{-36,-10},{36,-32}}, textString=
                                              "O2"),
        Line(
          points={{0,-50},{0,-90}},
          thickness=0.5),
        Line(points={{50,0},{88,0}}),
        Text(extent={{-80,100},{80,60}}, textString=
                                             "%name")}));
end sensor_O2;

model sensor_Q
    "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput Q
    annotation (Placement(transformation(
        origin={0,-98},
        extent={{-10,-10},{10,10}},
        rotation=270)));
equation

  In.Q + Out.Q = 0;
  Q = In.Q;
  // eventually abs(In.Q) to be shure to have pos. signal
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
"));
end sensor_Q;

model sensor_TKN "Ideal TKN and total nitrogen sensor"

  extends WasteWater.Icons.sensor_TKN;
  extends Interfaces.stoichiometry;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput TKN[2]
    annotation (Placement(transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0.0;
  TKN[1] = In.Snh + In.Snd + In.Xnd + i_xb*(In.Xbh + In.Xba) + i_xp*(In.Xp + In.Xi);
  TKN[2] = TKN[1] + In.Sno;

  annotation (
    Documentation(info="This component measures the Total Kjeldal Nitrogen (TKN) and the
total nitrogen (N_total) concentration [g/m3] of ASM1 wastewater
and provides the result as output signal (to be further processed
with blocks of the Modelica.Blocks library).

signal[1] - TKN
signal[2] - N_total
"));
end sensor_TKN;

model sensor_TSS
    "Ideal sensor to measure total suspended solids concentration (ASM1)"

  extends WasteWater.Icons.sensor_TSS;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput ME annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0;

  TSS = 0.75*(In.Xs + In.Xbh + In.Xba + In.Xp + In.Xi);
  // the factor 0.75 needs to be adapted due to plant dependency
  // 0.75 is from the COST Benchmark configuration

  annotation (
    Documentation(info="This component measures the total suspended solids concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_TSS;

  model sensor_AE "Areation Energy"

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
"),   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
  end sensor_AE;


  model sensor_SP "Sludge Production"

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

    der(SP*T) = 1/1000 * 0.75 * (In.Xs + (u.Xi +In.Xi) + In.Xbh + In.Xba)*In.Q;

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
"),   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics));
  end sensor_SP;

  model sensor_EQ "Effluent Quality"

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
    Modelica.Blocks.Interfaces.RealOutput EQ(start=0) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={0,-100})));
    Real T(start=1e-3);
  equation

    In.Q + Out.Q = 0;
    // eventually abs(In.Q) to be shure to have pos. signal

    S_Nkj = In.Snh + In.Snd + In.Xnd + i_xb*(In.Xbh + In.Xba) + i_xp*(In.Xp + In.Xi);
    SS_e = 0.75*(In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp);
    BOD = 0.25*(In.Ss + In.Xs + (1-f_p)*(In.Xbh + In.Xba));
    COD = In.Ss + In.Si + In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp;
    der(T) = 1.0;
    der(EQ*T) =1/1000*(2 * SS_e + 1 * COD + 30 * S_Nkj + 10 * In.Sno + 2 * BOD)  * In.Q;
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
"),   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics));
  end sensor_EQ;

  model sensor_IQ "Influent Quality"

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
            extent={{-10,-10},{10,10}},
          rotation=270,
          origin={0,-100})));
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
"),   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics));
  end sensor_IQ;

  model sensor_PE "Pumping Energy"

    extends WasteWater.Icons.sensor_Q;
    Modelica.Blocks.Interfaces.RealInput Qa annotation (Placement(transformation(extent={{-110,
              -50},{-90,-30}})));
    Modelica.Blocks.Interfaces.RealInput Qr annotation (Placement(transformation(extent={{-110,
              -10},{-90,10}})));
    Modelica.Blocks.Interfaces.RealInput Qw annotation (Placement(transformation(extent={{-110,50},
              {-90,30}})));
    Modelica.Blocks.Interfaces.RealOutput PE(start=0) annotation (Placement(transformation(extent={{92,-10},
              {112,10}})));
    Real T(start=1e-3);
  equation
    der(T) = 1.0;
    der(PE*T) = 0.004 * Qa + 0.008 * Qr + 0.05 * Qw;

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end sensor_PE;

  model sensor_ME "Me"

    extends WasteWater.Icons.sensor_TSS;
    Modelica.Blocks.Interfaces.RealInput Kla1 annotation (Placement(transformation(extent={{-110,
              100},{-90,80}})));
    Modelica.Blocks.Interfaces.RealInput Kla2 annotation (Placement(transformation(extent={{-110,62},
              {-90,42}})));
    Modelica.Blocks.Interfaces.RealInput Kla3 annotation (Placement(transformation(extent={{-110,10},
              {-90,-10}})));
    Modelica.Blocks.Interfaces.RealInput Kla4 annotation (Placement(transformation(extent={{-110,
              -40},{-90,-60}})));
    Modelica.Blocks.Interfaces.RealInput Kla5 annotation (Placement(transformation(extent={{-110,-80},
              {-90,-100}})));

    Modelica.Blocks.Interfaces.RealOutput ME(start=0) annotation (Placement(
          transformation(extent={{88,-10},{108,10}})));

    Real From_ax_1;
    Real From_ax_2;
    Real From_ax_3;
    Real From_ax_4;
    Real From_ax_5;
    Real T(start=1e-3);
    Real V = 1333;

  equation
    der(T) = 1;
    From_ax_1 = if (Kla3 < 20) then 0.005*V else 0;
    From_ax_2 = if (Kla4 < 20) then 0.005*V else 0;
    From_ax_3 = if (Kla5 < 20) then 0.005*V else 0;
    From_ax_4 = if (Kla5 < 20) then 0.005*V else 0;
    From_ax_5 = if (Kla5 < 20) then 0.005*V else 0;

    der(ME*T) = 24*(From_ax_1 + From_ax_2 + From_ax_3 + From_ax_4 + From_ax_5);

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end sensor_ME;

  model sensor_OCI "Overall cost index"

    extends WasteWater.Icons.sensor_Q;
    /*Modelica.Blocks.Interfaces.RealInput EC annotation (Placement(transformation(extent={{-110,
            -36},{-90,-56}})));*/
    Modelica.Blocks.Interfaces.RealInput AE annotation (Placement(transformation(extent={{-110,84},
              {-90,104}})));
    Modelica.Blocks.Interfaces.RealInput PE annotation (Placement(transformation(extent={{-110,
              -104},{-90,-84}})));
    Modelica.Blocks.Interfaces.RealInput SP annotation (Placement(transformation(extent={{-110,
              -30},{-90,-50}})));
    Modelica.Blocks.Interfaces.RealInput ME annotation (Placement(transformation(extent={{-110,50},
              {-90,30}})));
    Modelica.Blocks.Interfaces.RealOutput OCI annotation (Placement(transformation(extent={{92,-10},
              {112,10}})));
  equation
    //OCI = AE + PE + 5*SP +3*EC + ME;
    OCI = AE + PE + 5*SP +3*0 + ME;
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end sensor_OCI;

  model sensor_EC "unused"

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
  end sensor_EC;

  package Examples
    "Demonstration examples of the components of the ASM1 library"

    extends Modelica.Icons.Library;

    class BenchPlant "COST Benchmark WWTP Configuration"
      import WasteWater;

      //Q_air=34574.2654508612 is equal to a Kla of 10 h^-1 from COST benchmark
      //Q_air=12100.99290780142 is equal to a Kla of 3.5 h^-1 from COST benchmark
      extends Modelica.Icons.Example;

      WasteWater.BSM1.EffluentSink Effluent
        annotation (Placement(transformation(extent={{217,28},{237,48}})));
      WasteWater.BSM1.SludgeSink WasteSludge
        annotation (Placement(transformation(extent={{279,-46},{299,-26}})));
      WasteWater.BSM1.SecClarModTakacs Settler
        annotation (Placement(transformation(extent={{159,24},{179,44}})));
      WasteWater.BSM1.nitri tank4(V=1333)
        annotation (Placement(transformation(extent={{11,23},{31,43}})));
      WasteWater.BSM1.nitri tank3(V=1333)
        annotation (Placement(transformation(extent={{-35,23},{-15,43}})));
      WasteWater.BSM1.deni tank2
        annotation (Placement(transformation(extent={{-84,23},{-64,43}})));
      WasteWater.BSM1.deni tank1
        annotation (Placement(transformation(extent={{-125,24},{-105,44}})));
      Modelica.Blocks.Sources.Constant Constant1(k=55338)
                                                 annotation (Placement(
            transformation(extent={{-45,-13},{-65,7}})));
      Modelica.Blocks.Sources.Constant Constant2(k=18446)
                                                 annotation (Placement(
            transformation(extent={{36,-43},{16,-23}})));
      WasteWater.ASM1.divider2 divider
        annotation (Placement(transformation(extent={{112,23},{132,43}})));
      WasteWater.ASM1.mixer3 mixer
        annotation (Placement(transformation(extent={{-164,25},{-144,45}})));

      WasteWater.BSM1.WWSource Source
        annotation (Placement(transformation(extent={{-237,78},{-217,98}})));
      WasteWater.BSM1.sensor_IQ IQ
        annotation (Placement(transformation(extent={{-207,91},{-187,71}})));
      WasteWater.BSM1.sensor_EQ EQ
        annotation (Placement(transformation(extent={{186,50},{206,30}})));
      WasteWater.BSM1.sensor_SP SP
        annotation (Placement(transformation(extent={{255,-47},{275,-27}})));
      WasteWater.BSM1.sensor_Q Qw
        annotation (Placement(transformation(extent={{223,-51},{243,-31}})));
      WasteWater.BSM1.sensor_Q Qr
        annotation (Placement(transformation(extent={{142,-23},{122,-3}})));
      WasteWater.BSM1.sensor_Q Qa
        annotation (Placement(transformation(extent={{128,6},{108,26}})));
      WasteWater.BSM1.sensor_PE PE
        annotation (Placement(transformation(extent={{238,-81},{258,-61}})));
      WasteWater.BSM1.sensor_ME ME
        annotation (Placement(transformation(extent={{83,45},{103,65}})));
      WasteWater.BSM1.sensor_AE AE
        annotation (Placement(transformation(extent={{83,66},{103,86}})));
      WasteWater.BSM1.sensor_OCI OCI
        annotation (Placement(transformation(extent={{317,-2},{337,18}})));
      WasteWater.BSM1.pump pumpA(Q_max=60000)
        annotation (Placement(transformation(extent={{-84,22},{-104,2}})));
      WasteWater.BSM1.pump pumpR
        annotation (Placement(transformation(extent={{2,-7},{-18,-27}})));
      Modelica.Blocks.Sources.Constant Constant3(k=385)
                                                 annotation (Placement(
            transformation(extent={{151,-53},{171,-33}})));
      WasteWater.BSM1.pump pump
        annotation (Placement(transformation(extent={{194,-28},{214,-48}})));
      WasteWater.BSM1.nitri_tank5 tank5
        annotation (Placement(transformation(extent={{83,19},{103,39}})));
      Modelica.Blocks.Sources.CombiTimeTable Input(
        tableOnFile=true,
        tableName="Inf_dry",
        columns=integer(({2,3,4,5,6,7,8,9,10,11,12,13,14,15})),
        fileName=ModelicaServices.ExternalReferences.loadResource(
            "modelica://WasteWater/Resources/ASM1/Inf_dry.txt"))
        annotation (Placement(transformation(extent={{-266,76},{-246,96}})));
    equation
      connect(tank3.Out, tank4.In) annotation (Line(points={{-15,33},{-15,33},{
              11,33}}));
      connect(tank3.In, tank2.Out) annotation (Line(points={{-35,33},{-35,33},{
              -64,33}}));
      connect(tank1.Out, tank2.In) annotation (Line(points={{-105,34},{-105,33},
              {-84,33}}));

      connect(divider.Out1, Settler.Feed) annotation (Line(
          points={{132,35.6},{144,35.6},{144,35.4},{159,35.4}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(mixer.Out, tank1.In) annotation (Line(
          points={{-144,34.6},{-129,34.6},{-129,34},{-125,34}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(Source.Out, IQ.In)
        annotation (Line(points={{-217.2,81},{-207,81}}, smooth=Smooth.None));
      connect(Settler.Effluent, EQ.In) annotation (Line(points={{179.2,39.7},{
              181.6,39.7},{181.6,40},{186,40}}, smooth=Smooth.None));
      connect(EQ.Out, Effluent.In)
        annotation (Line(points={{206,40},{217,40}}, smooth=Smooth.None));
      connect(SP.Out, WasteSludge.In) annotation (Line(
          points={{275,-37},{270,-37},{270,-37.2},{279,-37.2}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(SP.u, Settler.Return) annotation (Line(
          points={{254.8,-31.8},{254.8,-33},{255,-33},{255,-13},{153,-13},{153,
              24.4},{166,24.4}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(Qw.Out, SP.In) annotation (Line(
          points={{243,-41},{249,-41},{249,-41.2},{255.2,-41.2}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(PE.Qw, Qw.Q) annotation (Line(
          points={{238,-67},{233,-67},{233,-50.8}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(PE.Qr, Qr.Q) annotation (Line(
          points={{238,-71},{132,-71},{132,-22.8}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(PE.Qa, Qa.Q) annotation (Line(
          points={{238,-75},{118,-75},{118,6.2}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(IQ.Out, mixer.In1) annotation (Line(
          points={{-187,81},{-176,81},{-176,38.5},{-164,38.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(ME.Kla4, tank4.Kla) annotation (Line(
          points={{83,50},{21.2,50},{21.2,36}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ME.Kla3, tank3.Kla) annotation (Line(
          points={{83,55},{-24.8,55},{-24.8,36}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ME.Kla2, tank2.Kla) annotation (Line(
          points={{83,60.2},{-74,60.2},{-74,38.4}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ME.Kla1, tank1.Kla) annotation (Line(
          points={{83,64},{-115,64},{-115,39.4}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(AE.Kla4, tank4.Kla) annotation (Line(
          points={{83,76},{21,76},{21,50},{21.2,50},{21.2,36}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(AE.Kla3, tank3.Kla) annotation (Line(
          points={{83,80},{-25,80},{-25,55},{-24.8,55},{-24.8,36}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(SP.SP, OCI.SP) annotation (Line(
          points={{265,-46.8},{265,-56},{301,-56},{301,4},{317,4}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(PE.PE, OCI.PE) annotation (Line(
          points={{258.2,-71},{310,-71},{310,-1.4},{317,-1.4}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ME.ME, OCI.ME) annotation (Line(
          points={{102.8,55},{301,55},{301,12},{317,12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(AE.AE, OCI.AE) annotation (Line(
          points={{103.2,76},{310,76},{310,17.4},{317,17.4}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Constant1.y, pumpA.u) annotation (Line(
          points={{-66,-3},{-74,-3},{-74,9.5},{-85.1,9.5}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(pumpR.u, Constant2.y) annotation (Line(
          points={{0.9,-19.5},{5,-19.5},{5,-20},{10,-20},{10,-33},{15,-33}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pumpA.Out, mixer.In3) annotation (Line(
          points={{-104,9.2},{-169,9.2},{-169,30.5},{-164,30.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(pumpR.Out, mixer.In2) annotation (Line(
          points={{-18,-19.8},{-183,-19.8},{-183,34.5},{-164,34.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(pump.Out, Qw.In) annotation (Line(
          points={{214,-40.8},{219,-40.8},{219,-41},{223,-41}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Settler.Waste, pump.In) annotation (Line(
          points={{172,24.4},{176,24.4},{176,24},{184,24},{184,-34.7},{194,
              -34.7}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Constant3.y, pump.u) annotation (Line(
          points={{172,-43},{184,-43},{184,-40.5},{195.1,-40.5}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(tank4.Out, tank5.In) annotation (Line(
          points={{31,33},{49,33},{49,29},{83,29}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(tank5.Out, divider.In) annotation (Line(
          points={{103,29},{93,29},{93,33.3},{112,33.3}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(ME.Kla5, tank5.Kla) annotation (Line(
          points={{83,46},{83,32},{93,32}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(AE.Kla5, tank5.Kla) annotation (Line(
          points={{83,72},{83,32},{93,32}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Input.y,Source. data) annotation (Line(
          points={{-245,86},{-248,86},{-248,88},{-236,88}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Qa.In, divider.Out2) annotation (Line(points={{128,16},{142,16},{
              142,31.5},{132,31.5}}, color={0,0,255}));
      connect(Qa.Out, pumpA.In) annotation (Line(points={{108,16},{92,16},{92,
              15.3},{-84,15.3}}, color={0,0,0}));
      connect(Qr.In, Settler.Return) annotation (Line(points={{142,-13},{145,
              -13},{145,-12},{152,-12},{152,-11},{153,-11},{153,24.4},{166,24.4}},
            color={0,0,255}));
      connect(Qr.Out, pumpR.In)
        annotation (Line(points={{122,-13},{2,-13},{2,-13.7}}, color={0,0,0}));
      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-270,-110},{380,110}},
            grid={1,1})),
        Documentation(info="This ASM1 plant consists of 2 denitrification tanks (tank1 and tank2),
3 nitrification tanks (tank3 - tank5) and a secondary clarifier by Takacs.
Furthermore there are 2 control loops modelled.
This configuration corresponds to the COST simulation benchmark [1].

Change into the directory ../ASM1 and translate the model.
Before simulating the model load initial values from the script file bench_asm1.mos
that is provided besides the model.
A 14 days dynamic influent data file is provided. So you may simulate up to 14 days.
But start with 1 day as it may take some time for simulation.
After simulation you may have a look at internal concentrations but most interesting
are the relevant concentrations at the effluent of a plant which can be viewed via the
sensors at the effluent of the secondary clarifier.

References:

[1] J.B. Copp: The COST Simulation Benchmark. 2000. http://www.ensic.u-nancy.fr/COSTWWTP/


PS: For those who want to reproduce the exact figures from the COST simulation benchmark some remarks:
    The aeration system in this library is different from that in COST, so be sure to produce an airflow
    corresponding to the desired Kla in COST. Furthermore in this library biological parameters are standard
    parameters from the ASM1 distribution and implemented with temperature dependency which may vary a bit from
    the parameter set used in COST.
    But it is possible. During the validation phase of this library the steady state and dynamic results
    from the COST simulation benchmark could exactly be reproduced.
"),     experiment(
          StopTime=14,
          __Dymola_NumberOfIntervals=100,
          __Dymola_Algorithm="lsodar"),
        __Dymola_experimentSetupOutput);
    end BenchPlant;

    model JenaSecClarModTakacs
      "Secondary Clarifier Model based on Takacs prepared for a special plant"

      extends WasteWater.Icons.SecClar;
      extends SecClar.Takacs.Interfaces.ratios;
      package SCP = SecClar.Takacs;
      import SI = Modelica.SIunits;
      package WI = Interfaces;
      package WWU = WasteWater.WasteWaterUnits;
      parameter SI.Length hsc=4.0 "height of secondary clarifier";
      parameter Integer n=10 "number of layers of SC model";
      parameter SI.Length zm=hsc/(1.0*n)
        "height of m-th secondary clarifier layer";
      parameter SI.Area Asc=1500.0 "area of secondary clarifier";
      parameter WWU.MassConcentration Xt=3000.0 "threshold for X";
      // total sludge concentration in clarifier feed
      WWU.MassConcentration Xf;
      // layers 1 to 10
      WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                {-90,24}})));
      WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                47},{112,67}})));
      WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                -106},{-20,-86}})));
      WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                -106},{40,-86}})));
      SCP.bottom_layer S1(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        rXs=rXs,
        rXbh=rXbh,
        rXba=rXba,
        rXp=rXp,
        rXi=rXi,
        rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
      SCP.lower_layer S2(
        zm=zm,
        Asc=Asc,
        Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
      SCP.lower_layer S3(
        zm=zm,
        Asc=Asc,
        Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
      SCP.feed_layer S4(
        zm=zm,
        Asc=Asc,
        Xf=Xf) annotation (Placement(transformation(extent={{-36,-36},{34,-22}})));
      SCP.upper_layer S5(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-36,-16},{34,-4}})));
      SCP.upper_layer S6(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-36,2},{34,16}})));
      SCP.upper_layer S7(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
      SCP.upper_layer S8(
        zm=zm,
        Asc=Asc,
        Xt=Xt,
        Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
      SCP.upper_layer S9(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
      SCP.top_layer S10(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt,
        rXs=rXs,
        rXbh=rXbh,
        rXba=rXba,
        rXp=rXp,
        rXi=rXi,
        rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
    equation

      connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
              -2.22045e-15,-74}}));
      connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
              -2.22045e-15,-55}}));
      connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
              -2.22045e-15,40}}));
      connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
              -2.22045e-15,78}}));
      connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
              -2.22045e-15,59}}));
      connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{30,-93},{30,-96}}));
      connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
              67.5,57},{102,57}}));
      connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-30,-93},{-30,
              -96}}));
      connect(S4.Dn, S3.Up) annotation (Line(points={{-1,-36},{-1,-38},{0,-38},{0,
              -40}}));
      connect(S4.Up, S5.Dn) annotation (Line(points={{-1,-22},{-2,-22},{-2,-22},{
              0,-22},{-1,-22},{-1,-16}}));
      connect(S5.Up, S6.Dn) annotation (Line(points={{-1,-4},{-2,-4},{-2,-4},{-2,
              -2},{-1,-2},{-1,2}}));
      connect(S6.Up, S7.Dn) annotation (Line(points={{-1,16},{-1,18},{0,18},{0,21}}));
      connect(Feed, S4.In) annotation (Line(points={{-100,14},{-67.5,14},{-67.5,
              -28.72},{-36,-28.72}}));

      // total sludge concentration in clarifier feed
      Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

      // ratios of solid components
      rXs = Feed.Xs/Xf;
      rXbh = Feed.Xbh/Xf;
      rXba = Feed.Xba/Xf;
      rXp = Feed.Xp/Xf;
      rXi = Feed.Xi/Xf;
      rXnd = Feed.Xnd/Xf;

      annotation (
        Documentation(info="This component models an ASM1 10 - layer secondary clarifier with 6 layers above the feed_layer (including top_layer)
and 3 layers below the feed_layer (including bottom_layer).

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  Xt  -  threshold value for Xtss [mg/l]
"));
    end JenaSecClarModTakacs;
    annotation (
      Documentation(info="This package contains example ASM1 wastewater treatment plant models to demonstrate the usage of
the WasteWater.ASM1 library.
Open the models and simulate them according to the description provided in the models.

The following demo models are present:

 - SmallPlant
 - BenchPlant
 - ComplexPlant

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
  end Examples;

  package Interfaces
    "Connectors and partial ASM1 model for Wastewater Treatment Modelling"

    extends Modelica.Icons.Library;

    connector WWFlowAsm1in "Inflow connector of ASM1 components"
      package WWU = WasteWater.WasteWaterUnits;
      flow WWU.VolumeFlowRate Q;
      WWU.MassConcentration Si;
      WWU.MassConcentration Ss;
      WWU.MassConcentration Xi;
      WWU.MassConcentration Xs;
      WWU.MassConcentration Xbh;
      WWU.MassConcentration Xba;
      WWU.MassConcentration Xp;
      WWU.MassConcentration So;
      WWU.MassConcentration Sno;
      WWU.MassConcentration Snh;
      WWU.MassConcentration Snd;
      WWU.MassConcentration Xnd;
      WWU.Alkalinity Salk;

      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={0,0,191},
              fillPattern=FillPattern.Solid), Text(
              extent={{-88,92},{88,-94}},
              lineColor={255,255,255},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString=
                   "%name")}),
        Documentation(info="Connectors WWFlowAsm1in and WWFlowAsm1out are nearly identical.
The difference is in the icons to more easily identify the inflow and outflow
side of a component.
The connector consists of one flow variable and 13 potential variables (ASM1 concentrations).
"));

    end WWFlowAsm1in;

    connector WWFlowAsm1out "Outflow connector of ASM1 components"
      package WWU = WasteWater.WasteWaterUnits;
      flow WWU.VolumeFlowRate Q;
      WWU.MassConcentration Si;
      WWU.MassConcentration Ss;
      WWU.MassConcentration Xi;
      WWU.MassConcentration Xs;
      WWU.MassConcentration Xbh;
      WWU.MassConcentration Xba;
      WWU.MassConcentration Xp;
      WWU.MassConcentration So;
      WWU.MassConcentration Sno;
      WWU.MassConcentration Snh;
      WWU.MassConcentration Snd;
      WWU.MassConcentration Xnd;
      WWU.Alkalinity Salk;
      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-88,92},{94,-92}}, textString=
                             "%name")}),
        Documentation(info="Connectors WWFlowAsm1in and WWFlowAsm1out are nearly identical.
The difference is in the icons to more easily identify the inflow and outflow
side of a component.
The connector consists of one flow variable and 13 potential variables (ASM1 concentrations).
"));
    end WWFlowAsm1out;

    connector AirFlow "Airflow connector"
      package WWU = WasteWater.WasteWaterUnits;

      flow WWU.VolumeFlowRate Q_air;
      annotation (
        Documentation(info="The Airflow connector consists of a flow variable describing the exchange of
air between blower and nitrification tank."));
    end AirFlow;

    partial model stoichiometry "ASM1 stoichiometric coefficients"
      // Stoichiometric parameters based on the original ASM1 publication//
      parameter Real Y_h=0.67
        "Heterotrophic Yield [g Xbh COD formed/(g COD utilised)]";
      parameter Real Y_a=0.24
        "Autotrophic Yield [g Xba COD formed/(g N utilised)]";
      parameter Real f_p=0.08 "Fraction of biomass to particulate products [-]";
      parameter Real i_xb=0.08 "Fraction nitrogen in biomass [g N/(g COD)]";
      parameter Real i_xp=0.06
        "Fraction nitrogen in particulate products [g N/(g COD)]";
      annotation (
        Documentation(info=
              "This is a partial model providing the stoichiometric coefficients of the ASM1 model."));
    end stoichiometry;

    partial model ASM1base "Base class of WWTP modelling by ASM1"
      extends Interfaces.stoichiometry;
      package WWU = WasteWater.WasteWaterUnits;

      // parameters based on the original ASM1 publication based on 15 deg C
      Real mu_h "Maximum heterotrophic growth rate f(T) [day^-1]";
      Real b_h "Heterotrophic decay rate f(T) [day^-1]";
      Real mu_a "Maximum autotrophic growth rate f(T) [day^-1]";
      //Real K_nh "Half-saturation (auto. growth) f(T) [g NH-N/m3]";
      Real b_a "Autotrophic decay rate f(T) [day^-1]";
      Real k_a "Ammonification rate f(T) [m3/(g COD day)]";
      Real k_h "Maximum specific hydrolysis rate f(T) [g Xs/(g Xbh COD day)]";
      Real K_x "Half-saturation (hydrolysis) f(T) [g Xs/(g Xbh COD)]";
      parameter Real mu_h_T=4.0
        "Maximum heterotrophic growth rate at T=15 deg C [day^-1]";
      parameter Real b_h_T=0.3
        "Heterotrophic decay rate at T=15 deg C [day^-1]";
      parameter Real mu_a_T=0.5
        "Maximum autotrophic growth rate at T=15 deg C[day^-1]";
      parameter Real b_a_T=0.05 "Autotrophic decay rate at T=15 deg C [day^-1]";
      parameter Real k_a_T=0.05
        "Ammonification rate at T=15 deg C [m3/(g COD day)]";
      parameter Real k_h_T=3.0
        "Maximum specific hydrolysis rate at T=15 deg C [g Xs/(g Xbh COD day)]";
      parameter Real K_x_T=0.1
        "Half-saturation (hydrolysis) at T=15 deg C [g Xs/(g Xbh COD)]";
      parameter Real K_nh=1.0 "Half-saturation (auto. growth) [g NH-N/m3]";
      parameter Real K_s=10.0 "Half-saturation (hetero. growth) [g COD/m3]";
      parameter Real K_oh=0.2 "Half-saturation (hetero. oxygen) [g O/m3]";
      parameter Real K_no=0.5 "Half-saturation (nitrate) [g NO-N/m3]";
      parameter Real K_oa=0.4 "Half-saturation (auto. oxygen) [g O/m3]";
      parameter Real ny_g=0.8 "Anoxic growth rate correction factor [-]";
      parameter Real ny_h=0.8 "Anoxic hydrolysis rate correction factor [-]";
      WWU.MassConcentration Si(fixed=true) "Soluble inert organic matter";
      WWU.MassConcentration Ss(fixed=true) "Readily biodegradable substrate";
      WWU.MassConcentration Xi(fixed=true) "Particulate inert organic matter";
      WWU.MassConcentration Xs(fixed=true) "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh(fixed=true) "Active heterotrophic biomass";
      WWU.MassConcentration Xba(fixed=true) "Active autotrophic biomass";
      WWU.MassConcentration Xp(fixed=true)
        "Particulate products from biomass decay";
      WWU.MassConcentration So(fixed=true) "Dissolved oxygen";
      WWU.MassConcentration Sno(fixed=true) "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh(fixed=true) "Ammonium nitrogen";
      WWU.MassConcentration Snd(fixed=true)
        "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd(fixed=true)
        "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk(fixed=true) "Alkalinity";
      Real p1;
      Real p2;
      Real p3;
      Real p4;
      Real p5;
      Real p6;
      Real p7;
      Real p8;
      Real r1;
      Real r2;
      Real r3;
      Real r4;
      Real r5;
      Real r6;
      Real r7;
      Real r8;
      Real r9;
      Real r10;
      Real r11;
      Real r12;
      Real r13;
      Real inputSi;
      Real inputSs;
      Real inputXi;
      Real inputXs;
      Real inputXbh;
      Real inputXba;
      Real inputXp;
      Real inputSo;
      Real inputSno;
      Real inputSnh;
      Real inputSnd;
      Real inputXnd;
      Real inputSalk;
      Real aeration;

      Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{
                -110,-10},{-90,10}})));
      Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{
                90,-10},{110,10}})));
      Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
              extent={{50,40},{60,50}})));
      Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
              extent={{-110,30},{-90,50}})));
    equation

      // Temperature dependent Kinetic parameters based on 15 deg C //
      // may be adapted to 10 or 20 deg C

      mu_h =mu_h_T;
      b_h =b_h_T;
      mu_a =mu_a_T;
      //K_nh=1.0*exp(0.069*(T.signal[1]-15));
      b_a =b_a_T;
      k_a =k_a_T;
      k_h =k_h_T;
      K_x =K_x_T;

      // Process Rates //

      p1 = mu_h*(Ss/(K_s + Ss))*(So/(K_oh + So))*Xbh;
      p2 = mu_h*(Ss/(K_s + Ss))*(K_oh/(K_oh + So))*(Sno/(K_no + Sno))*ny_g*Xbh;
      p3 = mu_a*(Snh/(K_nh + Snh))*(So/(K_oa + So))*Xba;
      p4 = b_h*Xbh;
      p5 = b_a*Xba;
      p6 = k_a*Snd*Xbh;
      p7 = k_h*((Xs/Xbh)/(K_x + (Xs/Xbh)))*((So/(K_oh + So)) + ny_h*(K_oh/(K_oh
         + So))*(Sno/(K_no + Sno)))*Xbh;
      p8 = p7*Xnd/Xs;

      // biochemical reactions

      r1 = 0;
      r2 = (-p1 - p2)/Y_h + p7;
      r3 = 0;
      r4 = (1 - f_p)*(p4 + p5) - p7;
      r5 = p1 + p2 - p4;
      r6 = p3 - p5;
      r7 = f_p*(p4 + p5);
      r8 = -((1 - Y_h)/Y_h)*p1 - ((4.57 - Y_a)/Y_a)*p3;
      r9 = -((1 - Y_h)/(2.86*Y_h))*p2 + p3/Y_a;
      r10 = -i_xb*(p1 + p2) - (i_xb + (1/Y_a))*p3 + p6;
      r11 = -p6 + p8;
      r12 = (i_xb - f_p*i_xp)*(p4 + p5) - p8;
      r13 = -i_xb/14*p1 + ((1 - Y_h)/(14*2.86*Y_h) - (i_xb/14))*p2 - ((i_xb/14)
         + 1/(7*Y_a))*p3 + p6/14;

      // derivatives

      der(Si) = inputSi + r1;
      der(Ss) = inputSs + r2;
      der(Xi) = inputXi + r3;
      der(Xs) = inputXs + r4;
      der(Xbh) = inputXbh + r5;
      der(Xba) = inputXba + r6;
      der(Xp) = inputXp + r7;
      der(So) = inputSo + r8 + aeration;
      der(Sno) = inputSno + r9;
      der(Snh) = inputSnh + r10;
      der(Snd) = inputSnd + r11;
      der(Xnd) = inputXnd + r12;
      der(Salk) = inputSalk + r13;

      // Outputs

      Out.Q + In.Q = 0;
      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = Xi;
      Out.Xs = Xs;
      Out.Xbh = Xbh;
      Out.Xba = Xba;
      Out.Xp = Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = Xnd;
      Out.Salk = Salk;

      MeasurePort.Si = Si;
      MeasurePort.Ss = Ss;
      MeasurePort.Xi = Xi;
      MeasurePort.Xs = Xs;
      MeasurePort.Xbh = Xbh;
      MeasurePort.Xba = Xba;
      MeasurePort.Xp = Xp;
      MeasurePort.So = So;
      MeasurePort.Sno = Sno;
      MeasurePort.Snh = Snh;
      MeasurePort.Snd = Snd;
      MeasurePort.Xnd = Xnd;
      MeasurePort.Salk = Salk;

      annotation (
        Documentation(info="This partial model provides connectors and equations that are needed in the biological
components (nitrification and denitrification tank) for ASM1 wastewater treatment plant models.
Parameters are coded according the ASM1 [1,2] standard distribution.
Changes to this parameters are subject to the modeller.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de


References:

[1] M. Henze and C.P.L. Grady Jr and W. Gujer and G.v.R. Marais and T. Matsuo:
    Activated Sludge Model No.1. IAWQ, 1987.
[2] M. Henze and W.Gujer and T. Mino and. M.v. Loosdrecht: Activated Sludge
    Models ASM1, ASM2, ASM2d, and ASM3. IWA Task Group on Mathematical Modelling
    for Design and Operation of Biological Wastewater Treatment, 2000.


This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2002, Gerald Reichl
"));
    end ASM1base;
    annotation (
      Documentation(info="This package contains connectors and interfaces (partial models) for
wastewater treatment components based on the Acticated Sludge Model No.1 (ASM1).

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2001, Gerald Reichl
"));
  end Interfaces;

  package PreClar "Primary clarifier modelling based on ASM1"

    extends Modelica.Icons.Library;

    model preclar1 "Dynamic ASM1 Primary Clarifier Model"
      import Modelica.Math.log;

      // dynamic primary clarifier tank, based on Otterpohl
      // to be used for feed forward calculation, e.g. influent data needed

      package WWU = WasteWaterUnits;
      extends WasteWater.Icons.preclar1;

      // tank specific parameters

      parameter Modelica.SIunits.Volume V=500
        "Volume of primary clarifier tank";
      Real hrt_h "hydraulic residence time in primary sedimentation tank [h]";

        // Real hrt_min "hydraulic residence time in primary sedimentation tank [min]";
      Real n_COD "efficiency of COD removal [%]";
      Real n_X "efficiency transformed to particulate fractions [%]";
      WWU.MassConcentration Si "Soluble inert organic matter";
      WWU.MassConcentration Ss "Readily biodegradable substrate";
      WWU.MassConcentration Xi "Particulate inert organic matter";
      WWU.MassConcentration Xs "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh "Active heterotrophic biomass";
      WWU.MassConcentration Xba "Active autotrophic biomass";
      WWU.MassConcentration Xp "Particulate products from biomass decay";
      WWU.MassConcentration So "Dissolved oxygen";
      WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh "Ammonium nitrogen";
      WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk "Alkalinity";
      Real CODin;
      Real CODout;
      Real XCODin;
      Real H;
      BSM1.Interfaces.WWFlowAsm1in In
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      BSM1.Interfaces.WWFlowAsm1out Out
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      BSM1.Interfaces.WWFlowAsm1out MeasurePort
        annotation (Placement(transformation(extent={{32,90},{42,100}})));
    equation

      // calculation of the hydraulic residence time
      hrt_h = V/In.Q*24;
      //hrt_min = V/In.Q * 24 * 60;

        // n_COD according Otterpohl and Freund 1992 "Dynamic Models for Clarifiers"
      n_COD = 2.7*(log(hrt_h*hrt_h) + 9)/100;
      // n_COD according Otterpohl 1995, Dissertation
      // n_COD = (1.45 + 6.15 * log(hrt_min))/100;

      XCODin = In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;
      // particulate COD in the influent
      CODin = In.Si + In.Ss + XCODin;
      // total COD in the influent

      CODout = Out.Si + Out.Ss + Out.Xi + Out.Xs + Out.Xbh + Out.Xba + Out.Xp;

      H = n_COD*CODin/XCODin;

      // n_X can not be greater than 1
      // therefore is this check
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;

      // in this case the model needs to be modified by a new n_COD
      // n_COD_? = (2.88*XCODin/CODin - 0.118) * n_COD;

      // volume dependent dilution term of each concentration

      der(Si) = (In.Si - Si)*In.Q/V;
      der(Ss) = (In.Ss - Ss)*In.Q/V;
      der(Xi) = (In.Xi - Xi)*In.Q/V;
      der(Xs) = (In.Xs - Xs)*In.Q/V;
      der(Xbh) = (In.Xbh - Xbh)*In.Q/V;
      der(Xba) = (In.Xba - Xba)*In.Q/V;
      der(Xp) = (In.Xp - Xp)*In.Q/V;
      der(So) = (In.So - So)*In.Q/V;
      der(Sno) = (In.Sno - Sno)*In.Q/V;
      der(Snh) = (In.Snh - Snh)*In.Q/V;
      der(Snd) = (In.Snd - Snd)*In.Q/V;
      der(Xnd) = (In.Xnd - Xnd)*In.Q/V;
      der(Salk) = (In.Salk - Salk)*In.Q/V;

      // Outputs
      // this is just a reduction of particulate substances; n_X*X is not stored
      // so the amount of primary sludge removed is not calculated
      Out.Q + In.Q = 0;

      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = (1 - n_X)*Xi;
      Out.Xs = (1 - n_X)*Xs;
      Out.Xbh = (1 - n_X)*Xbh;
      Out.Xba = (1 - n_X)*Xba;
      Out.Xp = (1 - n_X)*Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = (1 - n_X)*Xnd;
      Out.Salk = Salk;

      MeasurePort.Si = Si;
      MeasurePort.Ss = Ss;
      MeasurePort.Xi = (1 - n_X)*Xi;
      MeasurePort.Xs = (1 - n_X)*Xs;
      MeasurePort.Xbh = (1 - n_X)*Xbh;
      MeasurePort.Xba = (1 - n_X)*Xba;
      MeasurePort.Xp = (1 - n_X)*Xp;
      MeasurePort.So = So;
      MeasurePort.Sno = Sno;
      MeasurePort.Snh = Snh;
      MeasurePort.Snd = Snd;
      MeasurePort.Xnd = (1 - n_X)*Xnd;
      MeasurePort.Salk = Salk;
      annotation (
        Documentation(info="This is an ASM1 dynamic primary clarifier model based on the theory
by Otterpohl and Freund.

Parameter:
  V - volume of the preclarifier
"));
    end preclar1;

    model preclar2 "Static ASM1 Primary Clarifier Model"
      // static primary clarifier tank, based on Otterpohl
      // to be used for feed forward calculation, e.g. influent data needed

      import Modelica.Math.log;

      package WWU = WasteWaterUnits;
      extends WasteWater.Icons.preclar2;

      // tank specific parameters

      parameter Modelica.SIunits.Volume V=500
        "Volume of primary clarifier tank";
      Real hrt_h "hydraulic residence time in primary sedimentation tank [h]";

        //Real hrt_min "hydraulic residence time in primary sedimentation tank [min]";
      Real n_COD "efficiency of COD removal [%]";
      Real n_X "efficiency transformed to particulate fractions [%]";

      WWU.MassConcentration Si "Soluble inert organic matter";
      WWU.MassConcentration Ss "Readily biodegradable substrate";
      WWU.MassConcentration Xi "Particulate inert organic matter";
      WWU.MassConcentration Xs "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh "Active heterotrophic biomass";
      WWU.MassConcentration Xba "Active autotrophic biomass";
      WWU.MassConcentration Xp "Particulate products from biomass decay";
      WWU.MassConcentration So "Dissolved oxygen";
      WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh "Ammonium nitrogen";
      WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk "Alkalinity";

      Real CODin;
      Real CODout;
      Real XCODin;
      Real H;
      BSM1.Interfaces.WWFlowAsm1in In
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      BSM1.Interfaces.WWFlowAsm1out Out
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      BSM1.Interfaces.WWFlowAsm1out MeasurePort
        annotation (Placement(transformation(extent={{32,90},{42,100}})));
    equation

      // calculation of the hydraulic residence time
      hrt_h = V/In.Q*24;
      //hrt_min = V/In.Q * 24 * 60;

        // n_COD according Otterpohl and Freund 1992 "Dynamic Models for Clarifiers"
      n_COD = 2.7*(log(hrt_h*hrt_h) + 9)/100;
      // n_COD according Otterpohl 1995, Dissertation
      // n_COD = (1.45 + 6.15 * log(hrt_min))/100;

      XCODin = In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;
      // particulate COD in the influent
      CODin = In.Si + In.Ss + XCODin;
      // total COD in the influent

      CODout = Out.Si + Out.Ss + Out.Xi + Out.Xs + Out.Xbh + Out.Xba + Out.Xp;

      H = n_COD*CODin/XCODin;
      // n_X can not be greater than 1
      // therefore is this check
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;
      // in this case the model needs to be modified by a new n_COD
      // n_COD_? = (2.88*XCODin/CODin - 0.118) * n_COD;

      // volume dependent dilution term of each concentration

      0 = (In.Si - Si)*In.Q/V;
      0 = (In.Ss - Ss)*In.Q/V;
      0 = (In.Xi - Xi)*In.Q/V;
      0 = (In.Xs - Xs)*In.Q/V;
      0 = (In.Xbh - Xbh)*In.Q/V;
      0 = (In.Xba - Xba)*In.Q/V;
      0 = (In.Xp - Xp)*In.Q/V;
      0 = (In.So - So)*In.Q/V;
      0 = (In.Sno - Sno)*In.Q/V;
      0 = (In.Snh - Snh)*In.Q/V;
      0 = (In.Snd - Snd)*In.Q/V;
      0 = (In.Xnd - Xnd)*In.Q/V;
      0 = (In.Salk - Salk)*In.Q/V;

      // Outputs
      // this is just a reduction of particulate substances; n_X*X is not stored
      // so the amount of primary sludge removed is not calculated
      Out.Q + In.Q = 0;
      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = (1 - n_X)*Xi;
      Out.Xs = (1 - n_X)*Xs;
      Out.Xbh = (1 - n_X)*Xbh;
      Out.Xba = (1 - n_X)*Xba;
      Out.Xp = (1 - n_X)*Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = (1 - n_X)*Xnd;
      Out.Salk = Salk;

      MeasurePort.Si = Si;
      MeasurePort.Ss = Ss;
      MeasurePort.Xi = (1 - n_X)*Xi;
      MeasurePort.Xs = (1 - n_X)*Xs;
      MeasurePort.Xbh = (1 - n_X)*Xbh;
      MeasurePort.Xba = (1 - n_X)*Xba;
      MeasurePort.Xp = (1 - n_X)*Xp;
      MeasurePort.So = So;
      MeasurePort.Sno = Sno;
      MeasurePort.Snh = Snh;
      MeasurePort.Snd = Snd;
      MeasurePort.Xnd = (1 - n_X)*Xnd;
      MeasurePort.Salk = Salk;

      annotation (
        Documentation(info="This is an ASM1 static primary clarifier model based on the theory
by Otterpohl and Freund.

Parameter:
  V - volume of the preclarifier

"));
    end preclar2;

    model preclar3 "Inverse ASM1 Static Primary Clarifier Model"
      // static primary clarifier tank

        // to be used for backward calculation, e.g. effluent concentration data needed
      // signals need to be in the secuence COD, Sno, Snh, pH in the inputtable

      import Modelica.Math.log;
      extends WasteWater.Icons.preclar2;

      package WWU = WasteWater.WasteWaterUnits;

        // Interfaces.MeasurePort MeasurePort annotation (extent=[32, 90; 42, 100]);

      // tank specific parameters
      parameter Modelica.SIunits.Volume V=500
        "Volume of primary clarifier tank";
      parameter Real aSi=5/100
        "Fraction of Si of the total COD in the influent";
      parameter Real aSs=15/100
        "Fraction of Ss of the total COD in the influent";
      parameter Real aXi=15/100
        "Fraction of Xi of the total COD in the influent";
      parameter Real aXs=45/100
        "Fraction of Xs of the total COD in the influent";
      parameter Real aXbh=20/100
        "Fraction of Xbh of the total COD in the influent";
      parameter Real aXba=0/100
        "Fraction of Xba of the total COD in the influent";
      parameter Real aXp=0/100
        "Fraction of Xp of the total COD in the influent";
      parameter Real aSo=0.0 "Dissolved oxygen in the inflow [mg/l]";
      parameter Real aSnd=1/100 "Fraction Snd of Ss in the influent";
      parameter Real aXnd=3/100 "Fraction Xnd of Xs in the influent";
      parameter Real n_corr=1.0 "Correction faktor for the efficiency function";

      Real hrt_h "hydraulic residence time in primary sedimentation tank [h]";

        //Real hrt_min "hydraulic residence time in primary sedimentation tank [min]";
      Real n_COD "efficiency of COD removal [%]";
      Real n_X "efficiency transformed to particulate fractions [%]";

      WWU.MassConcentration Si "Soluble inert organic matter";
      WWU.MassConcentration Ss "Readily biodegradable substrate";
      WWU.MassConcentration Xi "Particulate inert organic matter";
      WWU.MassConcentration Xs "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh "Active heterotrophic biomass";
      WWU.MassConcentration Xba "Active autotrophic biomass";
      WWU.MassConcentration Xp "Particulate products from biomass decay";
      WWU.MassConcentration So "Dissolved oxygen";
      WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh "Ammonium nitrogen";
      WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk "Alkalinity";
      Real COD;
      Real CODin;
      Real CODout;
      Real XCOD;
      Real H;
      BSM1.Interfaces.WWFlowAsm1in In
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      BSM1.Interfaces.WWFlowAsm1out Out
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      Modelica.Blocks.Interfaces.RealInput MeasurePort[4]
        annotation (Placement(transformation(
            origin={38,90},
            extent={{-10,-10},{10,10}},
            rotation=270)));
    equation

      // calculation of the hydraulic residence time
      hrt_h = V/In.Q*24;
      //hrt_min = V/In.Q * 24 * 60;

        // n_COD according Otterpohl and Freund 1992 "Dynamic Models for Clarifiers"
      n_COD = n_corr*2.7*(log(hrt_h*hrt_h) + 9)/100;
      // n_COD according Otterpohl 1995, Dissertation
      // n_COD = (1.45 + 6.15 * log(hrt_min))/100;

      XCOD = In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;
      // particulate COD in the influent
      COD = In.Si + In.Ss + XCOD;
      // total COD in the influent

      CODin =MeasurePort[1]/(1 - n_COD);
      // total COD in the influent
      // above two CODs sould be the same

      CODout = Out.Si + Out.Ss + Out.Xi + Out.Xs + Out.Xbh + Out.Xba + Out.Xp;
      // this should be the same as MeasurePort.signal[1]

      H = n_COD*COD/XCOD;
      // n_X can not be greater than 1
      // therefor this check
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;
      // in this case the model needs to be modified by a new n_COD
      // n_COD_? = (2.88*XCODin/CODin - 0.118) * n_COD;

      // volume dependent dilution term of each concentration

      0 = (In.Si - Si)*In.Q/V;
      0 = (In.Ss - Ss)*In.Q/V;
      0 = (In.Xi - Xi)*In.Q/V;
      0 = (In.Xs - Xs)*In.Q/V;
      0 = (In.Xbh - Xbh)*In.Q/V;
      0 = (In.Xba - Xba)*In.Q/V;
      0 = (In.Xp - Xp)*In.Q/V;
      0 = (In.So - So)*In.Q/V;
      0 = (In.Sno - Sno)*In.Q/V;
      0 = (In.Snh - Snh)*In.Q/V;
      0 = (In.Snd - Snd)*In.Q/V;
      0 = (In.Xnd - Xnd)*In.Q/V;
      0 = (In.Salk - Salk)*In.Q/V;

      Out.Q + In.Q = 0;

      // Inputs
      In.Si = aSi*CODin;
      In.Ss = aSs*CODin;
      In.Xi = aXi*CODin;
      In.Xs = aXs*CODin;
      In.Xbh = aXbh*CODin;
      In.Xba = aXba*CODin;
      In.Xp = aXp*CODin;
      In.So = aSo;
      In.Sno =MeasurePort[2];
      In.Snh =MeasurePort[3];
      In.Snd = aSnd*In.Ss;
      In.Xnd = aXnd*In.Xs;
      In.Salk =1.8*exp(MeasurePort[4] - 6.4);

      // Outputs
      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = (1 - n_X)*Xi;
      Out.Xs = (1 - n_X)*Xs;
      Out.Xbh = (1 - n_X)*Xbh;
      Out.Xba = (1 - n_X)*Xba;
      Out.Xp = (1 - n_X)*Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = (1 - n_X)*Xnd;
      Out.Salk = Salk;

      annotation (
        Documentation(info="This is a special case of the ASM1 static primary clarifier model.
Here measurement data at the end (effluent) of the preclarifier needs to be provided.
This is typical for some real plants. Influent is then calculated.

Parameters:
  V   - volume of the preclarifier
  aS* - fractions of e.g. COD in influent (soluble components)
  aX* - fractions of e.g. COD in influent (particular components)
  n_corr- correction factor for the efficiency function

Dimension of InPort is 4.
  1 - Chemical Oxygen Demand (COD) at effluent of primary clarifier
  2 - nitrate nitrogen (Sno) at effluent of primary clarifier
  3 - ammonium nitrogen (Snh) at effluent of primary clarifier
  4 - pH-value at effluent of primary clarifier

"));
    end preclar3;
    annotation (
      Documentation(info="This package provides one dynamic and two static ASM1 primary clarifier
models based on Otterpohl [1].

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de


Reference:

[1] R. Otterpohl and M. Freund: Dynamic models for clarifier of activated sludge
    plants with dry and wet weather flows. Water Science and Technology. 26 (1992), pp. 1391-1400.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2001, Gerald Reichl
"));
  end PreClar;

  package SecClar "Library for secondary settling tank modelling based on ASM1"
  extends Modelica.Icons.Library;

    package Haertel "Secondary settling tank modelling by Haertel (ASM1)"
      extends Modelica.Icons.Library;

      package Interfaces
        "Connectors and partial models for ASM1 Secondary Clarifier Model by Haertel"

        extends Modelica.Icons.Library;

        connector UpperLayerPin "Connector above influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // effluent flow
          flow WWU.VolumeFlowRate Qe;

          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;

          // total sludge concentration (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers above the influent layer (feed_layer)."));
        end UpperLayerPin;

        connector LowerLayerPin "Connector below influent layer"
          package WWU = WasteWater.WasteWaterUnits;

          // return and waste sludge flow Qr, Qw
          flow WWU.VolumeFlowRate Qr;
          flow WWU.VolumeFlowRate Qw;
          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;
          // total sludge concentration in m-th layer
          WWU.MassConcentration X;

            // total sludge concentration and sink velocity in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;
          WWU.SedimentationVelocity vS_dn;
          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers below the influent layer (feed_layer)."));
        end LowerLayerPin;

        partial model SCParam "partial model providing clarifier parameters"
          import SI = Modelica.SIunits;
          package WWU = WasteWater.WasteWaterUnits;
          parameter SI.Length zm;
          parameter SI.Area Asc;
          parameter WWU.SludgeVolumeIndex ISV;
          annotation (
            Documentation(info="partial model providing clarifier parameters"));
        end SCParam;

        partial model SCVar "partial models providing variables"
          package WWU = WasteWater.WasteWaterUnits;
          WWU.MassConcentration X "total sludge concentration in m-th layer";
          WWU.SedimentationVelocity vS "sink velocity in m-th layer";
          WWU.SedimentationFlux Jsm "sedimentation flux m-th layer";

          WWU.MassConcentration Si "Soluble inert organic matter";
          WWU.MassConcentration Ss "Readily biodegradable substrate";
          WWU.MassConcentration So "Dissolved oxygen";
          WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
          WWU.MassConcentration Snh "Ammonium nitrogen";
          WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
          WWU.Alkalinity Salk "Alkalinity";
          annotation (
            Documentation(info="partial models providing variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"
          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;

        function vSfun "Sedimentation velocity function"

          // total sludge concentration in m-th layer in g/m3 or mg/l
          input Real X;
          //Sludge Volume Index
          input Real ISV;

          // sink velocity in m/d
          output Real vS;
        protected
          Real v0 "maximum settling velocity";
          Real nv "exponent as part of the Vesilind equation";
        algorithm
          v0 := (17.4*(exp(-0.0113*ISV)) + 3.931)*24;
          //[m/d]
          nv := (-0.9834*(exp(-0.00581*ISV)) + 1.043);
          //[l/g]
          vS := v0*exp(-nv*X/1000);
          annotation (
            Documentation(info="Sedimentation velocity function"));
        end vSfun;

        function omega "Omega correction function by Haertel"
          //vertical coordinate, bottom: z=0
          input Real z;
          // total sludge concentration in clarifier feed
          input Real Xf;
          //height of secondary clarifier
          input Real hsc;
          //height of m-th secondary clarifier layer
          input Real zm;
          //Sludge Volume Index
          input Real ISV;
          //number of layers above feed layer
          input Integer i;
          // correction function omega by Haertel based on [g/l]
          output Real omega;

        protected
          Real Xc "solids concentration at compression point";
          Real nv "exponent as part of the Vesilind equation";
          Real ht "height of transition point";
          Real hc "height of compressing point";
          Real B3;
          Real B4;

        algorithm
          Xc := 480/ISV;
          nv := 1.043 - 0.9834*exp(-0.00581*ISV);
          hc := (Xf/1000)*(hsc - zm*(i + 0.5))/Xc*(1.0 - 1.0/(Xc*nv));
          // unit change
          ht := min(2.0*hc, hsc - zm*(i + 0.5));

          B4 := 1.0 + 2.0*ISV/(100.0 + ISV);
          B3 := -((2*ISV + 100.0)/ISV)*hc^B4;

          omega := (1.0 - B3*ht^(-B4))/(1.0 - B3*z^(-B4));
          omega := min(1.0, omega);
          annotation (
            Documentation(info=
                  "This is Haertels omega correction function for the settling process."));
        end omega;
        annotation (
          Documentation(info="This package contains connectors and interfaces (partial models) for
the ASM1 secondary clarifier model based on Haertel [1] (settling velocity and omega correction function).

References:

[1] L. Haertel: Modellansaetze zur dynamischen Simulation des Belebtschlammverfahrens.
    TH Darmstadt, Dissertation, 1990.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2002, Gerald Reichl
"));
      end Interfaces;

      model SecClarModHaertel
        "ASM1 Secondary Settling Tank Model based on Haertel"
        import WasteWater;

        extends WasteWater.Icons.SecClar;
        extends Interfaces.ratios;
        package SCP = Haertel;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;
        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter Integer n=10 "number of layers of SC model";
        parameter SI.Length zm=hsc/(1.0*n)
          "height of m-th secondary clarifier layer";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.SludgeVolumeIndex ISV=130 "Sludge Volume Index";
        parameter Integer i=2
          "number of layers above current feed layer in this model";

        // total sludge concentration in clarifier feed
        WWU.MassConcentration Xf;

        // layers 1 to 10
        SCP.bottom_layer S1(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
        SCP.lower_layer S2(
          hsc=hsc,
          zm=zm,
          z=(zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
        SCP.lower_layer S3(
          hsc=hsc,
          zm=zm,
          z=(2*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
        SCP.lower_layer S4(
          hsc=hsc,
          zm=zm,
          z=(3*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
        SCP.lower_layer S5(
          hsc=hsc,
          zm=zm,
          z=(4*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
        SCP.lower_layer S6(
          hsc=hsc,
          zm=zm,
          z=(5*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
        SCP.lower_layer S7(
          hsc=hsc,
          zm=zm,
          z=(6*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
        SCP.feed_layer S8(
          hsc=hsc,
          zm=zm,
          z=(7*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
        SCP.upper_layer S9(
          zm=zm,
          Asc=Asc,
          ISV=ISV) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
        SCP.top_layer S10(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
                -2.22045e-15,-74}}));
        connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
                -2.22045e-15,-55}}));
        connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
                -2.22045e-15,-36}}));
        connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
                -2.22045e-15,2}}));
        connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
                -2.22045e-15,21}}));
        connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
                -2.22045e-15,40}}));
        connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
                -2.22045e-15,78}}));
        connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
                -2.22045e-15,-17}}));
        connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
                -2.22045e-15,59}}));
        connect(Feed, S8.In) annotation (Line(points={{-98,14},{-98,47.8},{-33,47.8}}));
        connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
                -100}}));
        connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
                67.5,57},{100,57}}));
        connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
                -100}}));

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

        annotation (
          Documentation(info="This component models an ASM1 10 - layer secondary clarifier with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Haertel`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  ISV -  Sludge Volume Index [ml/g]
  i   -  number of layers above feed layer
"));
      end SecClarModHaertel;

      model bottom_layer "Bottom layer of Haertel`s SC model"
        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;

        BSM1.Interfaces.WWFlowAsm1out PQr
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));
        BSM1.Interfaces.WWFlowAsm1out PQw
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);

        // sedimentation flux in bottom layer
        Jsm = 0.0;

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // upward connection
        Up.vS_dn = vS;
        Up.X_dn = X;

        // return and waste sludge volume flow rates
        PQr.Q + Up.Qr = 0;
        PQw.Q + Up.Qw = 0;

        // return sludge flow, solid and soluble components (ASM1)
        PQr.Si = Si;
        PQr.Ss = Ss;
        PQr.Xi = rXi*X;
        PQr.Xs = rXs*X;
        PQr.Xbh = rXbh*X;
        PQr.Xba = rXba*X;
        PQr.Xp = rXp*X;
        PQr.So = So;
        PQr.Sno = Sno;
        PQr.Snh = Snh;
        PQr.Snd = Snd;
        PQr.Xnd = rXnd*X;
        PQr.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        PQw.Si = Si;
        PQw.Ss = Ss;
        PQw.Xi = rXi*X;
        PQw.Xs = rXs*X;
        PQw.Xbh = rXbh*X;
        PQw.Xba = rXba*X;
        PQw.Xp = rXp*X;
        PQw.So = So;
        PQw.Sno = Sno;
        PQw.Snh = Snh;
        PQw.Snd = Snd;
        PQw.Xnd = rXnd*X;
        PQw.Salk = Salk;

        annotation (
          Documentation(info="This class models the lowest layer of an ASM1 secondary clarifier based on Haertel.

No sedimentation flux (mass exchange) with underneath but hydraulic and sedimentation flux (same direction)
with above layer.
From here return and waste sludge is removed.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end bottom_layer;

      model lower_layer "Layer below influent of Haertel`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";

        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;

        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if vS < Dn.vS_dn then omega*(vS*X) else omega*min(vS*X, Dn.vS_dn*Dn.
          X_dn);

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // downward connections
        Dn.Qr + Up.Qr = 0;
        Dn.Qw + Up.Qw = 0;

        Dn.X = X;
        Dn.SedFlux = -Jsm;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        // upward connections
        Up.vS_dn = vS;
        Up.X_dn = X;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the lowest layer (bottom_layer)
of an ASM1 secondary clarifier based on Haertel.

Hydraulic and sedimentation flux (mass exchange in same direction) with above and underneath layer.

Sedimentation flux is calculated based on the sedimentation velocity
function and the omega correction function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end lower_layer;

      model feed_layer "Influent layer of Haertel`s SC model"
        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;

        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";

        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;

        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        BSM1.Interfaces.WWFlowAsm1in In
          annotation (Placement(transformation(extent={{-110,-6},{-90,14}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if vS < Dn.vS_dn then omega*(vS*X) else omega*min(vS*X, Dn.vS_dn*Dn.
          X_dn);

        // ODE of solid component
        der(X) = (In.Q/Asc*Xf - (-Up.Qe)/Asc*X - (-(Dn.Qr + Dn.Qw))/Asc*X + Up.
          SedFlux - Jsm)/zm;

        // ODE of soluble components
        der(Si) = (In.Q*In.Si - (-Up.Qe)*Si - (-(Dn.Qr + Dn.Qw))*Si)/(Asc*zm);
        der(Ss) = (In.Q*In.Ss - (-Up.Qe)*Ss - (-(Dn.Qr + Dn.Qw))*Ss)/(Asc*zm);
        der(So) = (In.Q*In.So - (-Up.Qe)*So - (-(Dn.Qr + Dn.Qw))*So)/(Asc*zm);
        der(Sno) = (In.Q*In.Sno - (-Up.Qe)*Sno - (-(Dn.Qr + Dn.Qw))*Sno)/(Asc*zm);
        der(Snh) = (In.Q*In.Snh - (-Up.Qe)*Snh - (-(Dn.Qr + Dn.Qw))*Snh)/(Asc*zm);
        der(Snd) = (In.Q*In.Snd - (-Up.Qe)*Snd - (-(Dn.Qr + Dn.Qw))*Snd)/(Asc*zm);
        der(Salk) = (In.Q*In.Salk - (-Up.Qe)*Salk - (-(Dn.Qr + Dn.Qw))*Salk)/(Asc*
          zm);

        // volume flow rates
        In.Q + Up.Qe + Dn.Qr + Dn.Qw = 0;

        Dn.SedFlux = -Jsm;
        Dn.X = X;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the influent layer of an ASM1 secondary clarifier based on Haertel.

It receives the wastewater stream from the biological part (feed).
Hydraulic and sedimentation flux (mass exchange in opposite directions) with above layer
and hydraulic and sedimentation flux (mass exchange in same direction) with underneath layer.

Sedimentation flux is calculated based on the sedimentation velocity
function and the omega correction function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end feed_layer;

      model upper_layer "Layer above influent of Haertels`s SC"

        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = vS*X;

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // upward connections
        Up.Qe + Dn.Qe = 0;

        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the effluent layer (top_layer)
of an ASM1 secondary clarifier based on Haertel.

Hydraulic and sedimentation flux (mass exchange in opposite directions) with above and underneath layer.

Sedimentation flux is calculated based on the sedimentation velocity
function by Haertel."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end upper_layer;

      model top_layer "Effluent layer of Haertel`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        BSM1.Interfaces.WWFlowAsm1out Out
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = vS*X;

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // effluent volume flow rate
        Out.Q + Dn.Qe = 0;

        // effluent, solid and soluble components (ASM1)
        Out.Si = Si;
        Out.Ss = Ss;
        Out.Xi = rXi*X;
        Out.Xs = rXs*X;
        Out.Xbh = rXbh*X;
        Out.Xba = rXba*X;
        Out.Xp = rXp*X;
        Out.So = So;
        Out.Sno = Sno;
        Out.Snh = Snh;
        Out.Snd = Snd;
        Out.Xnd = rXnd*X;
        Out.Salk = Salk;

        annotation (
          Documentation(info="This class models the top layer of an ASM1 secondary clarifier based on Haertel.

No sedimentation flux (mass exchange) with above but hydraulic and sedimentation flux
(opposite directions) underneath.
From here effluent goes to the receiving water.

Sedimentation flux is calculated based on the sedimentation velocity
function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end top_layer;
      annotation (
        Documentation(info="This package contains classes (layer models) to built ASM1 secondary clarifier models, an Interfaces sub-library
and provides an ASM1 10-layer secondary clarifier model all bases on Haertel`s [1]
sedimentation velocity and omega correction functions.

A secondary clarifier layer model needs at least a top_layer, a feed_layer and a bottom_layer
and may have several upper_layer in between above the feed_layer and several lower_layer in
between below the feed_layer.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] L. Haertel: Modellansaetze zur dynamischen Simulation des Belebtschlammverfahrens.
    TH Darmstadt, Dissertation, 1990.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2002 - 2003, Gerald Reichl
"));
    end Haertel;

    package Krebs "Secondary settling tank modelling by Krebs (ASM1)"
      extends Modelica.Icons.Library;

      package Interfaces
        "Partial models for Secondary Clarifier Model by Krebs"

        extends Modelica.Icons.Library;

        partial model SCVar "partial models providing variables"
          package WWU = WasteWater.WasteWaterUnits;
          WWU.MassConcentration Xf "total sludge concentration";
          WWU.MassConcentration XB "sludge concentration in sludge layer";
          WWU.MassConcentration XR "sludge concentration of return";

          WWU.MassConcentration Si1(fixed=true)
            "Soluble inert organic matter in first stirrer tank of the excess layer";
          WWU.MassConcentration Ss1(fixed=true)
            "Readily biodegradable substrate in first stirrer tank of the excess layer";
          WWU.MassConcentration So1(fixed=true)
            "Dissolved oxygen in first stirrer tank of the excess layer";
          WWU.MassConcentration Sno1(fixed=true)
            "Nitrate and nitrite nitrogen in first stirrer tank of the excess layer";
          WWU.MassConcentration Snh1(fixed=true)
            "Ammonium nitrogen in first stirrer tank of the excess layer";
          WWU.MassConcentration Snd1(fixed=true)
            "Soluble biodegradable organic nitrogen in first stirrer tank of the excess layer";
          WWU.Alkalinity Salk1(fixed=true)
            "Alkalinity in first stirrer tank of the excess layer";

          WWU.MassConcentration Si2(fixed=true)
            "Soluble inert organic matter in second stirrer tank of the excess layer";
          WWU.MassConcentration Ss2(fixed=true)
            "Readily biodegradable substrate in second stirrer tank of the excess layer";
          WWU.MassConcentration So2(fixed=true)
            "Dissolved oxygen in second stirrer tank of the excess layer";
          WWU.MassConcentration Sno2(fixed=true)
            "Nitrate and nitrite nitrogen in second stirrer tank of the excess layer";
          WWU.MassConcentration Snh2(fixed=true)
            "Ammonium nitrogen in second stirrer tank of the excess layer";
          WWU.MassConcentration Snd2(fixed=true)
            "Soluble biodegradable organic nitrogen in second stirrer tank of the excess layer";
          WWU.Alkalinity Salk2(fixed=true)
            "Alkalinity in second stirrer tank of the excess layer";
          annotation (
            Documentation(info="partial models providing ASM1 variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"
          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;
        annotation (
          Documentation(info="This package contains partial models for ASM1 secondary clarifier models.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
      end Interfaces;

      model SecClarModKrebs "ASM1 Secondary Settling Tank Model based on Krebs"

        extends WasteWater.Icons.SecClarKrebs;
        import WWSC = WasteWater.BSM1.SecClar.Krebs.Interfaces;
        extends WWSC.SCVar;
        extends WWSC.ratios;

        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;
        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.SludgeVolumeIndex ISV=130 "Sludge Volume Index";
        Real te "thickening time in sludge layer in [d]";
        SI.Length hs "height of sludge layer";
        SI.Length he "height of excess layer";
        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

          //following expression is only for steady state initial equation of XB and is necessary

          //to calculate hs, if there would be problems with "initial()" in your modelica version
        //leave out this term and initial XB (or hs) manually e.g. via script-file
        if initial() then
          XB = Feed.Q/(0.7*(-(Return.Q + Waste.Q)))*Xf;
        end if;

        //thickening time in sludge layer in [d]
        te = 5/7*Asc*hs/(-(Return.Q + Waste.Q));

        //sludge concentration in sludge layer (unit of time in [h]) in [g/m3]
        XB = (1000/ISV*((te*24)^(1/3)))*1000;

        //sludge concentration of return
        XR = 0.7*XB;

        //ODE of height of sludge layer
        der(hs) = (Feed.Q*Xf - (-(Return.Q + Waste.Q))*XR)/(Asc/2*XB);

        //height of excess layer
        he = hsc - hs;

        // ODE of soluble components in first stirrer tank of the excess layer
        der(Si1) = (Feed.Q*Feed.Si - (-Effluent.Q)*Si1 - (-(Waste.Q + Return.Q))*
          Si1)/(Asc*he/2);
        der(Ss1) = (Feed.Q*Feed.Ss - (-Effluent.Q)*Ss1 - (-(Waste.Q + Return.Q))*
          Ss1)/(Asc*he/2);
        der(So1) = (Feed.Q*Feed.So - (-Effluent.Q)*So1 - (-(Waste.Q + Return.Q))*
          So1)/(Asc*he/2);
        der(Sno1) = (Feed.Q*Feed.Sno - (-Effluent.Q)*Sno1 - (-(Waste.Q + Return.Q))
          *Sno1)/(Asc*he/2);
        der(Snh1) = (Feed.Q*Feed.Snh - (-Effluent.Q)*Snh1 - (-(Waste.Q + Return.Q))
          *Snh1)/(Asc*he/2);
        der(Snd1) = (Feed.Q*Feed.Snd - (-Effluent.Q)*Snd1 - (-(Waste.Q + Return.Q))
          *Snd1)/(Asc*he/2);
        der(Salk1) = (Feed.Q*Feed.Salk - (-Effluent.Q)*Salk1 - (-(Waste.Q + Return.
          Q))*Salk1)/(Asc*he/2);

        // ODE of soluble components in second stirrer tank of the excess layer
        der(Si2) = ((-Effluent.Q)*Si1 - (-Effluent.Q)*Si2)/(Asc*he/2);
        der(Ss2) = ((-Effluent.Q)*Ss1 - (-Effluent.Q)*Ss2)/(Asc*he/2);
        der(So2) = ((-Effluent.Q)*So1 - (-Effluent.Q)*So2)/(Asc*he/2);
        der(Sno2) = ((-Effluent.Q)*Sno1 - (-Effluent.Q)*Sno2)/(Asc*he/2);
        der(Snh2) = ((-Effluent.Q)*Snh1 - (-Effluent.Q)*Snh2)/(Asc*he/2);
        der(Snd2) = ((-Effluent.Q)*Snd1 - (-Effluent.Q)*Snd2)/(Asc*he/2);
        der(Salk2) = ((-Effluent.Q)*Salk1 - (-Effluent.Q)*Salk2)/(Asc*he/2);

        // volume flow rates
        Feed.Q + Effluent.Q + Return.Q + Waste.Q = 0;

        // effluent, solid and soluble components (ASM1)
        Effluent.Si = Si2;
        Effluent.Ss = Ss2;
        Effluent.So = So2;
        Effluent.Sno = Sno2;
        Effluent.Snh = Snh2;
        Effluent.Snd = Snd2;
        Effluent.Salk = Salk2;
        Effluent.Xi = 0.0*XR;
        Effluent.Xs = 0.0*XR;
        Effluent.Xbh = 0.0*XR;
        Effluent.Xba = 0.0*XR;
        Effluent.Xp = 0.0*XR;
        Effluent.Xnd = 0.0*XR;

        // return sludge flow, solid and soluble components (ASM1)
        Return.Si = Si1;
        Return.Ss = Ss1;
        Return.So = So1;
        Return.Sno = Sno1;
        Return.Snh = Snh1;
        Return.Snd = Snd1;
        Return.Salk = Salk1;
        Return.Xi = rXi*XR;
        Return.Xs = rXs*XR;
        Return.Xbh = rXbh*XR;
        Return.Xba = rXba*XR;
        Return.Xp = rXp*XR;
        Return.Xnd = rXnd*XR;

        // waste sludge flow, solid and soluble components (ASM1)
        Waste.Si = Si1;
        Waste.Ss = Ss1;
        Waste.So = So1;
        Waste.Sno = Sno1;
        Waste.Snh = Snh1;
        Waste.Snd = Snd1;
        Waste.Salk = Salk1;
        Waste.Xi = rXi*XR;
        Waste.Xs = rXs*XR;
        Waste.Xbh = rXbh*XR;
        Waste.Xba = rXba*XR;
        Waste.Xp = rXp*XR;
        Waste.Xnd = rXnd*XR;
        annotation (
          Documentation(info="This component models an ASM1 secondary clarifier based on Krebs conceptional model.
It consists of two compartments: a \"sludge-bed\" and a clear water zone above.

Parameters:
  hsc -  height of clarifier [m]
  Asc -  surface area of secondary clarifier [m2]
  ISV -  Sludge Volume Index [ml/g]
"),       Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(
                extent={{-90,80},{92,14}},
                lineColor={0,0,255},
                lineThickness=0.5),
              Rectangle(
                extent={{-90,14},{92,-86}},
                lineColor={0,0,255},
                lineThickness=0.5),
              Polygon(
                points={{-8,-20},{-8,-38},{-16,-38},{0,-48},{16,-38},{8,-38},{8,-20},
                    {-8,-20}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,34},{-8,54},{-16,54},{0,64},{16,54},{8,54},{8,34},{-8,
                    34}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Text(extent={{-90,78},{-34,66}}, textString=
                                                   "top_layer"),
              Text(extent={{-90,20},{-30,-16}}, textString=
                                                    "bottom_layer"),
              Line(
                points={{-90,48},{92,48}},
                color={0,0,255},
                pattern=LinePattern.Dash)}));
      end SecClarModKrebs;
      annotation (
        Documentation(info="This package contains an ASM1 secondary clarifier model and an Interfaces sub-library
based on Krebs conceptional model [1].
The settler model consists of two compartments, a \"sludge-bed\" and a clear water zone above.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] P. Krebs and M. Armbruster and W. Rodi: Numerische Nachklaerbeckenmodelle. Korrespondenz Abwasser. 47 (7)
    2000. pp 985-999.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
    end Krebs;

    package Otterpohl "Secondary settling tank modelling by Otterpohl"
      extends Modelica.Icons.Library;

      package Interfaces
        "Connectors and partial ASM1 models for Secondary Clarifier Model by Otterpohl"

        extends Modelica.Icons.Library;

        connector UpperLayerPin "Connector above influent layer"

          package WWU = WasteWater.WasteWaterUnits;
          // effluent flow
          flow WWU.VolumeFlowRate Qe;
          // sedimentation flux (from micro and macro flocs)
          flow WWU.SedimentationFlux SedFlux_F;
          // caused by macro flocs
          flow WWU.SedimentationFlux SedFlux_S;
          // caused by micro flocs

            // sludge concentration of macro and micro flocs in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn_F;
          WWU.MassConcentration X_dn_S;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers above the influent layer (feed_layer)."));
        end UpperLayerPin;

        connector LowerLayerPin "Connector below influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // return and waste sludge flow Qr, Qw
          flow WWU.VolumeFlowRate Qr;
          flow WWU.VolumeFlowRate Qw;

          // sedimentation flux (from micro and macro flocs)
          flow WWU.SedimentationFlux SedFlux_F;
          // caused by macro flocs
          flow WWU.SedimentationFlux SedFlux_S;
          // caused by micro flocs

          // total sludge concentration of micro and macro flocs in m-th layer
          WWU.MassConcentration X_F;
          WWU.MassConcentration X_S;

            // total sludge concentration of micro and macro flocs in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn_F;
          WWU.MassConcentration X_dn_S;
          // sink velocity of macro flocs in (m-1)-th layer
          WWU.SedimentationVelocity vS_dn_F;
          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers below the influent layer (feed_layer)."));
        end LowerLayerPin;

        partial model SCParam "partial model providing clarifier parameters"
          import SI = Modelica.SIunits;
          package WWU = WasteWater.WasteWaterUnits;
          parameter SI.Length zm;
          parameter SI.Area Asc;
          parameter WWU.SludgeVolumeIndex ISV;
          parameter WWU.SedimentationVelocity vS_S=0.24;
          // 0.01[m/h]*24 -> [m/d]

          annotation (
            Documentation(info="partial model providing clarifier parameters"));
        end SCParam;

        partial model SCVar "partial models providing variables"

          package WWU = WasteWater.WasteWaterUnits;
          WWU.MassConcentration X "total sludge concentration in m-th layer";
          WWU.MassConcentration X_F "sludge concentration of macro flocs";
          WWU.MassConcentration X_S "sludge concentration of micro flocs";
          WWU.SedimentationVelocity vS_F "sink velocity of makro flocs";
          WWU.SedimentationFlux Jsm_F "sedimentation flux of macro flocs";
          WWU.SedimentationFlux Jsm_S "sedimentation flux of micro flocs";

          WWU.MassConcentration Si "Soluble inert organic matter";
          WWU.MassConcentration Ss "Readily biodegradable substrate";
          WWU.MassConcentration So "Dissolved oxygen";
          WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
          WWU.MassConcentration Snh "Ammonium nitrogen";
          WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
          WWU.Alkalinity Salk "Alkalinity";
          annotation (
            Documentation(info="partial models providing ASM1 variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"

          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;

        function vSfun "Sedimentation velocity function"

          // total sludge concentration in m-th layer in g/m3 or mg/l
          input Real X;
          //Sludge Volume Index
          input Real ISV;
          // sink velocity in m/d
          output Real vS;
        protected
          Real v0 "maximum settling velocity";
          Real nv "exponent as part of the Vesilind equation";
        algorithm
          v0 := (17.4*(exp(-0.0113*ISV)) + 3.931)*24;
          //[m/d]
          nv := (-0.9834*(exp(-0.00581*ISV)) + 1.043);
          //[l/g]
          vS := v0*exp(-nv*X/1000);
          annotation (
            Documentation(info="Sedimentation velocity function"));
        end vSfun;

        function omega "Omega correction function by Haertel"

          input Real z;
          //vertical coordinate, bottom: z=0
          input Real Xf;
          // total sludge concentration in clarifier feed
          input Real hsc;
          //height of secondary clarifier
          input Real zm;
          //height of m-th secondary clarifier layer
          input Real ISV;
          //Sludge Volume Index
          input Integer i;
          //number of layers above feed layer

          // correction function omega by Haertel based on [g/l]
          output Real omega;

        protected
          Real Xc "solids concentration at compression point";
          Real nv "exponent as part of the Vesilind equation";
          Real ht "height of transition point";
          Real hc "height of compressing point";
          Real B3;
          Real B4;

        algorithm
          Xc := 480/ISV;
          nv := 1.043 - 0.9834*exp(-0.00581*ISV);
          hc := (Xf/1000)*(hsc - zm*(i + 0.5))/Xc*(1.0 - 1.0/(Xc*nv));
          // unit change
          ht := min(2.0*hc, hsc - zm*(i + 0.5));

          B4 := 1.0 + 2.0*ISV/(100.0 + ISV);
          B3 := -((2*ISV + 100.0)/ISV)*hc^B4;

          omega := (1.0 - B3*ht^(-B4))/(1.0 - B3*z^(-B4));
          omega := min(1.0, omega);
          annotation (
            Documentation(info=
                  "This is Haertels omega correction function for the settling process."));
        end omega;
        annotation (
          Documentation(info="This package contains connectors and interfaces (partial models) for
the ASM1 secondary clarifier model based on Otterpohl [1] (two settling velocities for
distinction between micro and macro flocs and omega correction function).

References:

[1] R. Otterpohl and M. Freund: Dynamic models for clarifiers of activated sludge plants
    with dry and wet weather flows. Water Science and Technology. 26 (1992), pp 1391-1400.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
      end Interfaces;

      model SecClarModOtter
        "Secondary Clarifier Model based on Otterpohl (ASM1)"
        import WasteWater;

        extends WasteWater.Icons.SecClar;
        extends Interfaces.ratios;
        package SCP = Otterpohl;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;
        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter Integer n=10 "number of layers of SC model";
        parameter SI.Length zm=hsc/(1.0*n)
          "height of m-th secondary clarifier layer";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.SludgeVolumeIndex ISV=130 "Sludge Volume Index";
        parameter Integer i=2
          "number of layers above current feed layer in this model";

        // total sludge concentration in clarifier feed
        WWU.MassConcentration Xf;

        // layers 1 to 10
        SCP.bottom_layer S1(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
        SCP.lower_layer S2(
          hsc=hsc,
          zm=zm,
          z=(zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
        SCP.lower_layer S3(
          hsc=hsc,
          zm=zm,
          z=(2*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
        SCP.lower_layer S4(
          hsc=hsc,
          zm=zm,
          z=(3*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
        SCP.lower_layer S5(
          hsc=hsc,
          zm=zm,
          z=(4*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
        SCP.lower_layer S6(
          hsc=hsc,
          zm=zm,
          z=(5*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
        SCP.lower_layer S7(
          hsc=hsc,
          zm=zm,
          z=(6*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
        SCP.feed_layer S8(
          hsc=hsc,
          zm=zm,
          z=(7*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
        SCP.upper_layer S9(
          zm=zm,
          Asc=Asc,
          ISV=ISV) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
        SCP.top_layer S10(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
                -2.22045e-15,-74}}));
        connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
                -2.22045e-15,-55}}));
        connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
                -2.22045e-15,-36}}));
        connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
                -2.22045e-15,2}}));
        connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
                -2.22045e-15,21}}));
        connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
                -2.22045e-15,40}}));
        connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
                -2.22045e-15,78}}));
        connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
                -2.22045e-15,-17}}));
        connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
                -2.22045e-15,59}}));
        connect(Feed, S8.In) annotation (Line(points={{-98,14},{-98,47.8},{-33,47.8}}));
        connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
                -100}}));
        connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
                67.5,57},{100,57}}));
        connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
                -100}}));

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

        annotation (
          Documentation(info="This component models an ASM1 10 - layer secondary clarifier model with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Otterpohl`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  ISV -  Sludge Volume Index [ml/g]
  i   -  number of layers above feed layer
"));
      end SecClarModOtter;

      model bottom_layer "Bottom layer of Otterpohls`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        BSM1.Interfaces.WWFlowAsm1out PQr
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));
        BSM1.Interfaces.WWFlowAsm1out PQw
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);

        // sedimentation flux in bottom layer
        Jsm_F = 0.0;
        Jsm_S = 0.0;

        // ODE of solid component
        der(X_F) = ((Up.Qr + Up.Qw)/Asc*(Up.X_F - X_F) + Up.SedFlux_F)/zm;
        der(X_S) = ((Up.Qr + Up.Qw)/Asc*(Up.X_S - X_S) + Up.SedFlux_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // upward connection
        Up.vS_dn_F = vS_F;
        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;

        // return and waste sludge volume flow rates
        PQr.Q + Up.Qr = 0;
        PQw.Q + Up.Qw = 0;

        // return sludge flow, solid and soluble components (ASM1)
        PQr.Si = Si;
        PQr.Ss = Ss;
        PQr.Xi = rXi*X;
        PQr.Xs = rXs*X;
        PQr.Xbh = rXbh*X;
        PQr.Xba = rXba*X;
        PQr.Xp = rXp*X;
        PQr.So = So;
        PQr.Sno = Sno;
        PQr.Snh = Snh;
        PQr.Snd = Snd;
        PQr.Xnd = rXnd*X;
        PQr.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        PQw.Si = Si;
        PQw.Ss = Ss;
        PQw.Xi = rXi*X;
        PQw.Xs = rXs*X;
        PQw.Xbh = rXbh*X;
        PQw.Xba = rXba*X;
        PQw.Xp = rXp*X;
        PQw.So = So;
        PQw.Sno = Sno;
        PQw.Snh = Snh;
        PQw.Snd = Snd;
        PQw.Xnd = rXnd*X;
        PQw.Salk = Salk;

        annotation (
          Documentation(info="This class models the lowest layer of an ASM1 secondary clarifier based on Otterpohl.

No sedimentation flux (mass exchange) with underneath but hydraulic and sedimentation flux (same direction)
with above layer.
From here return and waste sludge is removed.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end bottom_layer;

      model lower_layer "Layer below influent of Otterpohl`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";
        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = if vS_F < Dn.vS_dn_F then omega*(vS_F*X_F) else omega*min(vS_F*X_F,
            Dn.vS_dn_F*Dn.X_dn_F);
        Jsm_S = omega*min(vS_S*X_S, vS_S*Dn.X_dn_S);

        // ODE of solid component
        der(X_F) = ((Up.Qr + Up.Qw)/Asc*(Up.X_F - X_F) + Up.SedFlux_F - Jsm_F)/zm;
        der(X_S) = ((Up.Qr + Up.Qw)/Asc*(Up.X_S - X_S) + Up.SedFlux_S - Jsm_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // downward connections
        Dn.Qr + Up.Qr = 0;
        Dn.Qw + Up.Qw = 0;

        Dn.X_F = X_F;
        Dn.X_S = X_S;
        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        // upward connections
        Up.vS_dn_F = vS_F;
        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the lowest layer (bottom_layer)
of an ASM1 secondary clarifier based on Otterpohl.

Hydraulic and sedimentation flux (mass exchange in same direction) with above and underneath layer.

Sedimentation flux is calculated based on two sedimentation velocities
(for macro and micro flocs) and the omega correction function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end lower_layer;

      model feed_layer "Influent layer of Otterpohl`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";
        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;
        Real fl;

        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        BSM1.Interfaces.WWFlowAsm1in In
          annotation (Placement(transformation(extent={{-110,-6},{-90,14}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);
        fl = (9.4/ISV)*exp(-1.07*Xf/1000);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = if vS_F < Dn.vS_dn_F then omega*(vS_F*X_F) else omega*min(vS_F*X_F,
            Dn.vS_dn_F*Dn.X_dn_F);
        Jsm_S = omega*min(vS_S*X_S, vS_S*Dn.X_dn_S);

        // ODE of solid component
        der(X_F) = (In.Q/Asc*Xf*(1 - fl) - (-Up.Qe)/Asc*X_F - (-(Dn.Qr + Dn.Qw))/
          Asc*X_F + Up.SedFlux_F - Jsm_F)/zm;
        der(X_S) = (In.Q/Asc*Xf*fl - (-Up.Qe)/Asc*X_S - (-(Dn.Qr + Dn.Qw))/Asc*X_S
           + Up.SedFlux_S - Jsm_S)/zm;

        X = X_F + X_S;

        // ODE of soluble components
        der(Si) = (In.Q*In.Si - (-Up.Qe)*Si - (-(Dn.Qr + Dn.Qw))*Si)/(Asc*zm);
        der(Ss) = (In.Q*In.Ss - (-Up.Qe)*Ss - (-(Dn.Qr + Dn.Qw))*Ss)/(Asc*zm);
        der(So) = (In.Q*In.So - (-Up.Qe)*So - (-(Dn.Qr + Dn.Qw))*So)/(Asc*zm);
        der(Sno) = (In.Q*In.Sno - (-Up.Qe)*Sno - (-(Dn.Qr + Dn.Qw))*Sno)/(Asc*zm);
        der(Snh) = (In.Q*In.Snh - (-Up.Qe)*Snh - (-(Dn.Qr + Dn.Qw))*Snh)/(Asc*zm);
        der(Snd) = (In.Q*In.Snd - (-Up.Qe)*Snd - (-(Dn.Qr + Dn.Qw))*Snd)/(Asc*zm);
        der(Salk) = (In.Q*In.Salk - (-Up.Qe)*Salk - (-(Dn.Qr + Dn.Qw))*Salk)/(Asc*
          zm);

        // volume flow rates
        In.Q + Up.Qe + Dn.Qr + Dn.Qw = 0;

        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;
        Dn.X_F = X_F;
        Dn.X_S = X_S;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the influent layer of an ASM1 secondary clarifier based on Otterpohl.

It receives the wastewater stream from the biological part (feed).
Hydraulic and sedimentation flux (mass exchange in opposite directions) with above layer
and hydraulic and sedimentation flux (mass exchange in same direction) with underneath layer.

Sedimentation flux is calculated based on two sedimentation velocities
(for macro and micro flocs) and the omega correction function by Haertel.
"),       Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end feed_layer;

      model upper_layer "Layer above influent of Otterpohl`s SC"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = vS_F*X_F;
        Jsm_S = vS_S*X_S;

        // ODE of solid component
        der(X_F) = (Dn.Qe/Asc*(Dn.X_dn_F - X_F) + Up.SedFlux_F - Jsm_F)/zm;
        der(X_S) = (Dn.Qe/Asc*(Dn.X_dn_S - X_S) + Up.SedFlux_S - Jsm_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;

        // upward connections
        Up.Qe + Dn.Qe = 0;

        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the effluent layer (top_layer)
of an ASM1 secondary clarifier based on Otterpohl.

Hydraulic and sedimentation flux (mass exchange in opposite directions) with above and underneath layer.

Sedimentation flux is calculated based on two sedimentation velocities
(for macro and micro flocs)."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end upper_layer;

      model top_layer "Effluent layer of Otterpohl`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        BSM1.Interfaces.WWFlowAsm1out Out
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = vS_F*X_F;
        Jsm_S = vS_S*X_S;

        // ODE of solid component
        der(X_F) = (Dn.Qe/Asc*(Dn.X_dn_F - X_F) - Jsm_F)/zm;
        der(X_S) = (Dn.Qe/Asc*(Dn.X_dn_S - X_S) - Jsm_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;

        // effluent volume flow rate
        Out.Q + Dn.Qe = 0;

        // effluent, solid and soluble components (ASM1)
        Out.Si = Si;
        Out.Ss = Ss;
        Out.Xi = rXi*X;
        Out.Xs = rXs*X;
        Out.Xbh = rXbh*X;
        Out.Xba = rXba*X;
        Out.Xp = rXp*X;
        Out.So = So;
        Out.Sno = Sno;
        Out.Snh = Snh;
        Out.Snd = Snd;
        Out.Xnd = rXnd*X;
        Out.Salk = Salk;
        annotation (
          Documentation(info="This class models the top layer of an ASM1 secondary clarifier based on Otterpohl.

No sedimentation flux (mass exchange) with above but hydraulic and sedimentation flux
(opposite directions) underneath.
From here effluent goes to the receiving water.

Sedimentation flux is calculated based on two sedimentation velocities
(for micro and macro flocs).
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end top_layer;
      annotation (
        Documentation(info="This package contains classes (layer models) to built ASM1 secondary clarifier models, an Interfaces sub-library
and provides an ASM1 10-layer secondary clarifier model all bases on Otterpohls`s [1]
sedimentation velocities for macro and micro flocs and the omega correction function.

A secondary clarifier layer model needs at least a top_layer, a feed_layer and a bottom_layer
and may have several upper_layer in between above the feed_layer and several lower_layer in
between below the feed_layer.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] R. Otterpohl and M. Freund: Dynamic models for clarifiers of activated sludge plants
    with dry and wet weather flows. Water Science and Technology. 26 (1992), pp 1391-1400.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
    end Otterpohl;

    package Simple "Simple ASM1 Secondary clarifier model"
      extends Modelica.Icons.Library;

      model SimpleSecClarMod "Simple ASM1 Secondary Clarifier Model"

        extends WasteWater.Icons.SecClarSimple;
        extends WasteWater.BSM1.SecClar.Takacs.Interfaces.ratios;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;

        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";

        WWU.MassConcentration Xf "total sludge concentration in clarifier feed";
        WWU.MassConcentration X "sludge concentration in clarifier";
        WWU.MassConcentration Si "Soluble inert organic matter";
        WWU.MassConcentration Ss "Readily biodegradable substrate";
        WWU.MassConcentration So "Dissolved oxygen";
        WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
        WWU.MassConcentration Snh "Ammonium nitrogen";
        WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
        WWU.Alkalinity Salk "Alkalinity";

        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

        // ODE of sludge concentration
        der(X) = (Feed.Q*Xf - (-(Waste.Q + Return.Q))*X)/(Asc*hsc);

        // ODE of soluble components
        der(Si) = (Feed.Q*Feed.Si - (-Effluent.Q)*Si - (-(Waste.Q + Return.Q))*Si)/
          (Asc*hsc);
        der(Ss) = (Feed.Q*Feed.Ss - (-Effluent.Q)*Ss - (-(Waste.Q + Return.Q))*Ss)/
          (Asc*hsc);
        der(So) = (Feed.Q*Feed.So - (-Effluent.Q)*So - (-(Waste.Q + Return.Q))*So)/
          (Asc*hsc);
        der(Sno) = (Feed.Q*Feed.Sno - (-Effluent.Q)*Sno - (-(Waste.Q + Return.Q))*
          Sno)/(Asc*hsc);
        der(Snh) = (Feed.Q*Feed.Snh - (-Effluent.Q)*Snh - (-(Waste.Q + Return.Q))*
          Snh)/(Asc*hsc);
        der(Snd) = (Feed.Q*Feed.Snd - (-Effluent.Q)*Snd - (-(Waste.Q + Return.Q))*
          Snd)/(Asc*hsc);
        der(Salk) = (Feed.Q*Feed.Salk - (-Effluent.Q)*Salk - (-(Waste.Q + Return.Q))
           *Salk)/(Asc*hsc);

        // volume flow rates
        Feed.Q + Effluent.Q + Return.Q + Waste.Q = 0;

        // effluent, solid and soluble components (ASM1)
        Effluent.Si = Si;
        Effluent.Ss = Ss;
        Effluent.Xi = 0.0*X;
        Effluent.Xs = 0.0*X;
        Effluent.Xbh = 0.0*X;
        Effluent.Xba = 0.0*X;
        Effluent.Xp = 0.0*X;
        Effluent.So = So;
        Effluent.Sno = Sno;
        Effluent.Snh = Snh;
        Effluent.Snd = Snd;
        Effluent.Xnd = 0.0*X;
        Effluent.Salk = Salk;

        // return sludge flow, solid and soluble components (ASM1)
        Return.Si = Si;
        Return.Ss = Ss;
        Return.Xi = rXi*X;
        Return.Xs = rXs*X;
        Return.Xbh = rXbh*X;
        Return.Xba = rXba*X;
        Return.Xp = rXp*X;
        Return.So = So;
        Return.Sno = Sno;
        Return.Snh = Snh;
        Return.Snd = Snd;
        Return.Xnd = rXnd*X;
        Return.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        Waste.Si = Si;
        Waste.Ss = Ss;
        Waste.Xi = rXi*X;
        Waste.Xs = rXs*X;
        Waste.Xbh = rXbh*X;
        Waste.Xba = rXba*X;
        Waste.Xp = rXp*X;
        Waste.So = So;
        Waste.Sno = Sno;
        Waste.Snh = Snh;
        Waste.Snd = Snd;
        Waste.Xnd = rXnd*X;
        Waste.Salk = Salk;

        annotation (
          Documentation(info="This component models very simple the secondary clarification process by
just using a single fully mixed tank which removes all particulate substances from the effluent
and returns the sludge. No sedimentation and compression, etc. is considered (for ASM1).

Parameters:
  hsc -    height of clarifier [m]
  Asc -    surface area of sec. clar. [m2]
"),       Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Polygon(
                points={{-20,-70},{20,-70},{4,-84},{-4,-84},{-20,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-4,-84},{4,-92}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-80,-48},{-36,-64},{38,-64},{80,-48},{-80,-48}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-80,62},{80,-40}},
                lineColor={223,191,159},
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Text(extent={{-80,98},{80,66}}, textString=
                                                  "%name"),
              Polygon(
                points={{-36,-64},{38,-64},{20,-70},{-20,-70},{-36,-64}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Line(
                points={{4,-92},{4,-84},{20,-70},{80,-48}},
                thickness=0.5),
              Rectangle(
                extent={{-80,-40},{80,-48}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{80,62},{92,54}},
                lineColor={0,127,255},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{80,54},{92,54}},
                thickness=0.5),
              Line(
                points={{-4,-92},{-4,-84},{-20,-70},{-80,-48},{-80,10}},
                thickness=0.5),
              Line(
                points={{-80,62},{-80,16}},
                thickness=0.5),
              Line(
                points={{-80,10},{-90,10}},
                thickness=0.5),
              Line(
                points={{-80,16},{-90,16}},
                thickness=0.5),
              Rectangle(
                extent={{-20,-92},{20,-98}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-20,-92},{-4,-92}},
                thickness=0.5),
              Line(
                points={{-20,-98},{20,-98}},
                thickness=0.5),
              Line(
                points={{20,-92},{4,-92}},
                thickness=0.5),
              Line(
                points={{80,-48},{80,54}},
                thickness=0.5),
              Text(extent={{-100,-60},{-40,-80}}, textString=
                                                      "return"),
              Text(extent={{40,-60},{100,-80}}, textString=
                                                    "waste"),
              Polygon(
                points={{16,44},{33,44},{31,52},{48,42},{31,31},{33,39},{16,39},{16,
                    44}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-46,32},{-29,32},{-31,40},{-14,30},{-31,19},{-29,27},{-46,
                    27},{-46,32}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{18,-26},{22,-26},{22,-42},{28,-40},{20,-54},{12,-40},{18,
                    -42},{18,-26}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={191,95,0},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-32,-10},{-28,-10},{-28,-26},{-22,-24},{-30,-38},{-38,-24},
                    {-32,-26},{-32,-10}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={191,95,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-90,16},{-80,10}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end SimpleSecClarMod;
      annotation (
        Documentation(info="This package just provides a very simple ASM1 secondary clarifier model
with no sludge storage, no sludge sedimentation and no use of layers.
The model consists of one tank removing all particulate substances.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2002, Gerald Reichl
"));
    end Simple;

    package Takacs "Secondary settling tank modelling by Takacs"
      extends Modelica.Icons.Library;

      package Interfaces
        "Connectors and partial models for Secondary Clarifier Model by Takacs"

        extends Modelica.Icons.Library;

        connector UpperLayerPin "Connector above influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // effluent flow
          flow WWU.VolumeFlowRate Qe;
          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;

            // total sludge concentration and sink velocity in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;
          WWU.SedimentationVelocity vS_dn;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers above the influent layer (feed_layer)."));
        end UpperLayerPin;

        connector LowerLayerPin "Connector below influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // return and waste sludge flow Qr, Qw
          flow WWU.VolumeFlowRate Qr;
          flow WWU.VolumeFlowRate Qw;

          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;

          // total sludge concentration in m-th layer
          WWU.MassConcentration X;

            // total sludge concentration and sink velocity in(m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;
          WWU.SedimentationVelocity vS_dn;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers below the influent layer (feed_layer)."));
        end LowerLayerPin;

        partial model SCParam "partial model providing clarifier parameters"
          import SI = Modelica.SIunits;
          parameter SI.Length zm;
          parameter SI.Area Asc;

          annotation (
            Documentation(info="partial model providing clarifier parameters"));
        end SCParam;

        partial model SCVar "partial models providing variables"

          package WWU = WasteWater.WasteWaterUnits;

          WWU.MassConcentration X "total sludge concentration in m-th layer";
          WWU.MassConcentration Xf
            "total sludge concentration in clarifier feed";
          WWU.SedimentationVelocity vS "sink velocity in m-th layer";
          WWU.SedimentationFlux Jsm "sedimentation flux m-th layer";

          WWU.MassConcentration Si "Soluble inert organic matter";
          WWU.MassConcentration Ss "Readily biodegradable substrate";
          WWU.MassConcentration So "Dissolved oxygen";
          WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
          WWU.MassConcentration Snh "Ammonium nitrogen";
          WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
          WWU.Alkalinity Salk "Alkalinity";
          annotation (
            Documentation(info="partial models providing ASM1 variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"

          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;

        function vSfun "Sedimentation velocity function"

          // total sludge concentration in m-th layer in g/m3 or mg/l
          input Real X;
          // total sludge concentration in clarifier feed in g/m3 or mg/l
          input Real Xf;
          // sink velocity in m/d
          output Real vS;

          parameter Real v0slash=250.0 "max. settling velocity in m/d";
          parameter Real v0=474.0 "max. Vesilind settl. veloc. in m/d";
          parameter Real rh=0.000576 "hindered zone settl. param. in m3/(g SS)";
          parameter Real rp=0.00286
            "flocculant zone settl. param. in m3/(g SS)";
          parameter Real fns=0.00228 "non-settleable fraction in -";
        algorithm

          // computation of sink velocity
          vS := max(0.0, min(v0slash, v0*(exp(-rh*(X - fns*Xf)) - exp(-rp*(X - fns*
            Xf)))));

          annotation (
            Documentation(info=
                  "Takacs double-exponential sedimentation velocity function."));
        end vSfun;
        annotation (
          Documentation(info="This package contains connectors and interfaces (partial models) for
the ASM1 secondary clarifier model based on Takacs [1] (double-exponential settling velocity).

References:

[1] I. Takacs and G.G. Patry and D. Nolasco: A dynamic model of the clarification-thickening
    process. Water Research. 25 (1991) 10, pp 1263-1271.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2001, Gerald Reichl
"));
      end Interfaces;

      model SecClarModTakacs "Secondary Clarifier ASM1 Model based on Takacs"

        extends WasteWater.Icons.SecClar;
        extends Interfaces.ratios;
        package SCP = Takacs;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;

        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter Integer n=10 "number of layers of SC model";
        parameter SI.Length zm=hsc/(1.0*n)
          "height of m-th secondary clarifier layer";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.MassConcentration Xt=3000.0 "threshold for X";

        // total sludge concentration in clarifier feed
        WWU.MassConcentration Xf;

        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));

        // layers 1 to 10
        SCP.bottom_layer S1(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
        SCP.lower_layer S2(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
        SCP.lower_layer S3(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
        SCP.lower_layer S4(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
        SCP.lower_layer S5(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
        SCP.feed_layer S6(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
        SCP.upper_layer S7(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          Xt=Xt) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
        SCP.upper_layer S8(
          zm=zm,
          Asc=Asc,
          Xt=Xt,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
        SCP.upper_layer S9(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          Xt=Xt) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
        SCP.top_layer S10(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          Xt=Xt,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
      equation

        connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
                -2.22045e-15,-74}}));
        connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
                -2.22045e-15,-55}}));
        connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
                -2.22045e-15,-36}}));
        connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
                -2.22045e-15,2}}));
        connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
                -2.22045e-15,21}}));
        connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
                -2.22045e-15,40}}));
        connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
                -2.22045e-15,78}}));
        connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
                -2.22045e-15,-17}}));
        connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
                -2.22045e-15,59}}));
        connect(Feed, S6.In) annotation (Line(points={{-100,10},{-67.5,10},{-67.5,
                9.8},{-35,9.8}}));
        connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
                -100}}));
        connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
                67.5,57},{100,57}}));
        connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
                -100}}));

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;
        annotation (
          Documentation(info="This component models an ASM1 10 - layer secondary clarifier model with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Takac`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  Xt  -  threshold value for Xtss [mg/l]

"));
      end SecClarModTakacs;

      model bottom_layer "Bottom layer of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        BSM1.Interfaces.WWFlowAsm1out PQr
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));
        BSM1.Interfaces.WWFlowAsm1out PQw
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);
        Jsm = 0.0;

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // upward connection
        Up.vS_dn = vS;
        Up.X_dn = X;

        // return and waste sludge volume flow rates
        PQr.Q + Up.Qr = 0;
        PQw.Q + Up.Qw = 0;

        // return sludge flow, solid and soluble components (ASM1)
        PQr.Si = Si;
        PQr.Ss = Ss;
        PQr.Xi = rXi*X;
        PQr.Xs = rXs*X;
        PQr.Xbh = rXbh*X;
        PQr.Xba = rXba*X;
        PQr.Xp = rXp*X;
        PQr.So = So;
        PQr.Sno = Sno;
        PQr.Snh = Snh;
        PQr.Snd = Snd;
        PQr.Xnd = rXnd*X;
        PQr.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        PQw.Si = Si;
        PQw.Ss = Ss;
        PQw.Xi = rXi*X;
        PQw.Xs = rXs*X;
        PQw.Xbh = rXbh*X;
        PQw.Xba = rXba*X;
        PQw.Xp = rXp*X;
        PQw.So = So;
        PQw.Sno = Sno;
        PQw.Snh = Snh;
        PQw.Snd = Snd;
        PQw.Xnd = rXnd*X;
        PQw.Salk = Salk;

        annotation (
          Documentation(info="This class models the lowest layer of an ASM1 secondary clarifier based on Takacs.

No sedimentation flux (mass exchange) with underneath but hydraulic and sedimentation flux (same direction)
with above layer.
From here return and waste sludge is removed.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end bottom_layer;

      model lower_layer "Layer below influent of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // downward connections
        Dn.Qr + Up.Qr = 0;
        Dn.Qw + Up.Qw = 0;

        Dn.X = X;
        Dn.SedFlux = -Jsm;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        // upward connections
        Up.vS_dn = vS;
        Up.X_dn = X;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the lowest layer (bottom_layer)
of an ASM1 secondary clarifier based on Takacs.

Hydraulic and sedimentation flux (mass exchange in same direction) with above and underneath layer.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end lower_layer;

      model feed_layer "Influent layer of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        BSM1.Interfaces.WWFlowAsm1in In
          annotation (Placement(transformation(extent={{-110,-6},{-90,14}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = (In.Q/Asc*Xf - (-Up.Qe)/Asc*X - (-(Dn.Qr + Dn.Qw))/Asc*X + Up.
          SedFlux - Jsm)/zm;

        // ODE of soluble components
        der(Si) = (In.Q*In.Si - (-Up.Qe)*Si - (-(Dn.Qr + Dn.Qw))*Si)/(Asc*zm);
        der(Ss) = (In.Q*In.Ss - (-Up.Qe)*Ss - (-(Dn.Qr + Dn.Qw))*Ss)/(Asc*zm);
        der(So) = (In.Q*In.So - (-Up.Qe)*So - (-(Dn.Qr + Dn.Qw))*So)/(Asc*zm);
        der(Sno) = (In.Q*In.Sno - (-Up.Qe)*Sno - (-(Dn.Qr + Dn.Qw))*Sno)/(Asc*zm);
        der(Snh) = (In.Q*In.Snh - (-Up.Qe)*Snh - (-(Dn.Qr + Dn.Qw))*Snh)/(Asc*zm);
        der(Snd) = (In.Q*In.Snd - (-Up.Qe)*Snd - (-(Dn.Qr + Dn.Qw))*Snd)/(Asc*zm);
        der(Salk) = (In.Q*In.Salk - (-Up.Qe)*Salk - (-(Dn.Qr + Dn.Qw))*Salk)/(Asc*
          zm);

        // volume flow rates
        In.Q + Up.Qe + Dn.Qr + Dn.Qw = 0;

        Dn.SedFlux = -Jsm;
        Dn.X = X;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        Up.vS_dn = vS;
        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;
        annotation (
          Documentation(info="This class models the influent layer of an ASM1 secondary clarifier based on Takacs.

It receives the wastewater stream from the biological part (feed).
Hydraulic and sedimentation flux (mass exchange in opposite directions) with above layer
and hydraulic and sedimentation flux (mass exchange in same direction) with underneath layer.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end feed_layer;

      model upper_layer "Layer above influent of Takac`s SC"
        // Xt = Xthreshold

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        parameter WWU.MassConcentration Xt;

        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if Dn.X_dn <= Xt then vS*X else min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // upward connections
        Up.Qe + Dn.Qe = 0;

        Up.vS_dn = vS;
        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the effluent layer (top_layer)
an ASM1 secondary clarifier based on Takacs.

Hydraulic and sedimentation flux (mass exchange in opposite directions) with above and underneath layer.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end upper_layer;

      model top_layer "Effluent layer of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;

        parameter WWU.MassConcentration Xt;
        // Xt = Xthreshold

        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        BSM1.Interfaces.WWFlowAsm1out Out
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if Dn.X_dn <= Xt then vS*X else min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // effluent volume flow rate
        Out.Q + Dn.Qe = 0;

        // effluent, solid and soluble components (ASM1)
        Out.Si = Si;
        Out.Ss = Ss;
        Out.Xi = rXi*X;
        Out.Xs = rXs*X;
        Out.Xbh = rXbh*X;
        Out.Xba = rXba*X;
        Out.Xp = rXp*X;
        Out.So = So;
        Out.Sno = Sno;
        Out.Snh = Snh;
        Out.Snd = Snd;
        Out.Xnd = rXnd*X;
        Out.Salk = Salk;
        annotation (
          Documentation(info="This class models the top layer of an ASM1 secondary clarifier based on Takacs.

No sedimentation flux (mass exchange) with above but hydraulic and sedimentation flux
(opposite directions) underneath.
From here effluent goes to the receiving water.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end top_layer;
      annotation (
        Documentation(info="This package contains classes (layer models) to built ASM1 secondary clarifier models, an Interfaces sub-library
and provides an ASM1 10-layer secondary clarifier model all bases on Takacs [1]
double exponential sedimentation velocity function.

A secondary clarifier layer model needs at least a top_layer, a feed_layer and a bottom_layer
and may have several upper_layer in between above the feed_layer and several lower_layer in
between below the feed_layer.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] I. Takacs and G.G. Patry and D. Nolasco: A dynamic model of the clarification-thickening
    process. Water Research. 25 (1991) 10, pp 1263-1271.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2002, Gerald Reichl
"));
    end Takacs;
  annotation (
    Documentation(info="This library provides a collection of ASM1 secondary clarifier models based on
several theories.

The library currently is structured in following sub-libraries.

 - Takacs: secondary clarifier model by Takacs et al [1]
 - Haertel: secondary clarifier model by Haertel [2]
 - Otterpohl: secondary clarifier model by Otterpohl [3]
 - Krebs: secondary clarifier model by Krebs [4]
 - Simple: very basic secondary clarifier model


Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] I. Takacs and G.G. Patry and D. Nolasco: A dynamic model of the clarification-thickening
    process. Water Research. 25 (1991) 10, pp 1263-1271.

[2] L. Haertel: Modellansaetze zur dynamischen Simulation des Belebtschlammverfahrens.
    TH Darmstadt, Dissertation, 1990.

[3] R. Otterpohl and M. Freund: Dynamic models for clarifiers of activated sludge plants
    with dry and wet weather flows. Water Science and Technology. 26 (1992), pp 1391-1400.

[4] P. Krebs and M. Armbruster and W. Rodi: Numerische Nachklaerbeckenmodelle. Korrespondenz Abwasser. 47 (7)
    2000. pp 985-999.


This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2003, Gerald Reichl
"));
  end SecClar;
annotation (
  Documentation(info="This library contains components to build models of biological municipal
wastewater treatment plants based on the Activated Sludge Model No.1 (ASM1) by the
International Association on Water Quality (IAWQ) [1,2].


The library currently is structured in following sub-libraries.

 - Interfaces (partial ASM1 models and connectors)
 - PreClar (primary clarifier models)
 - SecClar (several secondary settling tank models)
 - Examples (wastewater treatment plant models)

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de


References:

[1] M. Henze and C.P.L. Grady Jr and W. Gujer and G.v.R. Marais and T. Matsuo:
    Activated Sludge Model No.1. IAWQ, 1987.
[2] M. Henze and W.Gujer and T. Mino and. M.v. Loosdrecht: Activated Sludge
    Models ASM1, ASM2, ASM2d, and ASM3. IWA Task Group on Mathematical Modelling
    for Design and Operation of Biological Wastewater Treatment, 2000.


This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2002, Gerald Reichl
"));
end BSM1;
