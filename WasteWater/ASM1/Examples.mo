within WasteWater.ASM1;
package Examples "Demonstration examples of the components of the ASM1 library"

  extends Modelica.Icons.Library;

  class SmallPlant "Small WWTP Configuration"
    import WasteWater;
    extends Modelica.Icons.Example;

    //Q_air=12100.99290780142 is equal to a Kla of 3.5 h^-1 from COST benchmark
    //Q_air=34574.2654508612 is equal to a Kla of 10 h^-1 from COST benchmark

    ASM1.EffluentSink Effluent annotation (Placement(transformation(extent={{88,
              -28},{108,-8}})));
    ASM1.SludgeSink WasteSludge annotation (Placement(transformation(extent={{
              87,-51},{107,-31}})));
    ASM1.divider2 divider annotation (Placement(transformation(extent={{20,-6},
              {40,14}})));
    ASM1.nitri tank3(V=1333) annotation (Placement(transformation(extent={{-6,
              -6},{14,14}})));
    ASM1.nitri tank2(V=1333) annotation (Placement(transformation(extent={{-34,
              -6},{-14,14}})));
    ASM1.deni tank1(V=3000) annotation (Placement(transformation(extent={{-65,
              -6},{-45,14}})));
    ASM1.mixer3 mixer annotation (Placement(transformation(extent={{-104,22},{
              -84,42}})));
    Modelica.Blocks.Sources.CombiTimeTable CombiTableTime(
      fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/Inf_dry.txt"),
      table=[0,0; 1,1],
      columns=integer(({16,3,4,5,6,7,8,9,10,11,12,13,14,15})),
      tableName="Inf_dry",
      tableOnFile=("Inf_dry") <> "NoName") annotation (Placement(transformation(
            extent={{-114,78},{-94,98}})));
    ASM1.WWSource WWSource annotation (Placement(transformation(extent={{-88,78},
              {-68,98}})));
    ASM1.blower blower1(Q_max=34574.2654508612) annotation (Placement(
          transformation(extent={{-33,-62},{-13,-42}})));
    ASM1.blower blower2(Q_max=34574.2654508612) annotation (Placement(
          transformation(extent={{-6,-62},{14,-42}})));
    ASM1.sensor_O2 sensor_O2 annotation (Placement(transformation(extent={{0,24},
              {20,44}})));
    Modelica.Blocks.Math.Feedback Feedback annotation (Placement(transformation(
            extent={{62,40},{82,60}})));
    Modelica.Blocks.Continuous.PI PI1(T=0.001, k=500, initType=Modelica.Blocks.Types.Init.InitialState)
      annotation (Placement(transformation(extent={{88,40},{108,60}})));
    Modelica.Blocks.Sources.Constant Constant1 annotation (Placement(
          transformation(extent={{-67,-87},{-47,-67}})));
    ASM1.pump RecyclePump(Q_max=46115) annotation (Placement(transformation(
          origin={-84,-12},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    ASM1.pump ReturnPump(Q_max=9223) annotation (Placement(transformation(
          origin={26,-26},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    ASM1.pump WastePump(Q_max=193) annotation (Placement(transformation(extent=
              {{59,-55},{79,-35}})));
    Modelica.Blocks.Sources.Constant Constant2 annotation (Placement(
          transformation(extent={{22,-68},{42,-48}})));
    Modelica.Blocks.Sources.Constant Temperature(k=15)
      annotation (Placement(transformation(extent={{-94,50},{-82,62}})));
    sensor_NH sensor_NH1 annotation (Placement(transformation(extent={{64,15},{
              80,31}})));
    WasteWater.ASM1.sensor_NO sensor_NO1 annotation (Placement(transformation(
            extent={{81,15},{97,31}})));
    WasteWater.ASM1.sensor_TKN sensor_TKN1 annotation (Placement(transformation(
            extent={{97,14},{113,30}})));
    WasteWater.ASM1.sensor_COD sensor_COD1 annotation (Placement(transformation(
            extent={{97,-5},{113,11}})));
    Modelica.Blocks.Sources.Step OxygenSetpoint(height=1.5)
      annotation (Placement(transformation(extent={{37,40},{57,60}})));
    WasteWater.ASM1.SecClar.Krebs.SecClarModKrebs Settler annotation (Placement(
          transformation(extent={{48,-4},{68,16}})));
    WasteWater.ASM1.sensor_TSS sensor_TSS1 annotation (Placement(transformation(
            extent={{32,14},{49,30}})));
  equation
    connect(tank3.Out, divider.In) annotation (Line(points={{14,4},{17,4},{17,
            4.3},{20,4.3}}));
    connect(mixer.Out, tank1.In) annotation (Line(points={{-84,31.6},{-77,31.6},
            {-77,4},{-65,4}}));
    connect(mixer.In1, WWSource.Out) annotation (Line(points={{-104,35.5},{-104,
            74},{-68,74},{-68,81},{-68.2,81}}));
    connect(CombiTableTime.y, WWSource.data)
      annotation (Line(points={{-93,88},{-87,88}}));
    connect(blower2.AirOut, tank3.AirIn) annotation (Line(points={{3,-42},{3,
            -24},{4,-24},{4,-5.8}}));
    connect(Feedback.y, PI1.u) annotation (Line(points={{81,50},{86,50}}));
    connect(PI1.y, blower2.u) annotation (Line(points={{109,50},{114,50},{114,
            -84},{18,-84},{18,-55},{13.8,-55},{13.8,-55}}));
    connect(divider.Out2, RecyclePump.In) annotation (Line(points={{40,2.5},{44,
            2.5},{44,-8.7},{-74,-8.7}}));
    connect(RecyclePump.Out, mixer.In3) annotation (Line(points={{-94,-14.8},{
            -104,-14.8},{-104,27.5}}));
    connect(ReturnPump.Out, mixer.In2) annotation (Line(points={{16,-28.8},{
            15.5,-28.8},{15.5,-30},{-112,-30},{-112,31.5},{-104,31.5}}));
    connect(sensor_O2.So, Feedback.u2)
      annotation (Line(points={{19.8,34},{72,34},{72,42}}));
    connect(Temperature.y, tank1.T)
      annotation (Line(points={{-81.4,56},{-71,56},{-71,8},{-65,8}}, color={0,0,
            255}));
    connect(Temperature.y, tank2.T)
      annotation (Line(points={{-81.4,56},{-39,56},{-39,8},{-34,8}}, color={0,0,
            255}));
    connect(Temperature.y, tank3.T) annotation (Line(points={{-81.4,56},{-39,56},
            {-39,14},{-5.9,14},{-5.9,8},{-6,8}},   color={0,0,255}));
    connect(OxygenSetpoint.y, Feedback.u1)
      annotation (Line(points={{58,50},{64,50}}, color={0,0,255}));
    connect(Constant1.y, blower1.u) annotation (Line(points={{-46,-77},{-7.2,
            -77},{-7.2,-55},{-13.2,-55}}, color={0,0,255}));
    connect(WastePump.Out, WasteSludge.In) annotation (Line(points={{79,-42.2},
            {81,-42.2},{81,-42},{83,-42},{83,-42},{87,-42}}));
    connect(WastePump.u, Constant2.y)
      annotation (Line(points={{60.1,-42.5},{46,-42.5},{46,-58},{43,-58}},
                                                                     color={0,0,
            255}));
    connect(tank2.Out, tank3.In) annotation (Line(points={{-14,4},{-6,4}}));
    connect(tank1.Out, tank2.In) annotation (Line(points={{-45,4},{-34,4}}));
    connect(blower1.AirOut, tank2.AirIn) annotation (Line(points={{-24,-42},{
            -24,-24},{-24,-5.8},{-24,-5.8}}));
    connect(Constant1.y, RecyclePump.u) annotation (Line(points={{-46,-77},{-39,
            -77},{-39,-14.5},{-75.1,-14.5}},
                                         color={0,0,255}));
    connect(Settler.Effluent, Effluent.In) annotation (Line(points={{68.2,11.7},
            {78,11.7},{78,-16},{88,-16}}));
    connect(Settler.Return, ReturnPump.In) annotation (Line(points={{55,-3.6},{
            55,-22.7},{36,-22.7}}));
    connect(WastePump.In, Settler.Waste) annotation (Line(points={{59,-48.3},{
            52,-48.3},{52,-31},{61,-31},{61,-3.6}}));
    connect(sensor_NH1.In, Settler.Effluent) annotation (Line(points={{72,15},{
            72,11.7},{68.2,11.7}}));
    connect(sensor_NO1.In, Settler.Effluent) annotation (Line(points={{89,15},{
            89,11.7},{68.2,11.7}}));
    connect(sensor_TKN1.In, Settler.Effluent) annotation (Line(points={{105,14},
            {105,11.7},{68.2,11.7}}));
    connect(sensor_COD1.In, Settler.Effluent) annotation (Line(points={{105,-5},
            {105,11.7},{68.2,11.7}}));
    connect(Constant2.y, ReturnPump.u) annotation (Line(points={{43,-58},{48,
            -58},{48,-28.5},{34.9,-28.5}}, color={0,0,255}));
    connect(divider.Out1, Settler.Feed) annotation (Line(points={{40,6.6},{44.5,
            6.6},{44.5,7.4},{48,7.4}}));
    connect(tank3.MeasurePort, sensor_O2.In) annotation (Line(points={{9.5,8.5},
            {9.5,16},{10,16},{10,25}}));
    connect(sensor_TSS1.In, divider.Out1) annotation (Line(points={{40.5,14},{
            40.5,10.6},{40,10.6},{40,6}}));

    annotation (
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-120,-100},{120,105}},
          grid={1,1}), graphics={Line(points={{-22,58},{-22,58}}, color={0,0,
                255})}),
      Documentation(info="This fictitious plant provides an ASM1 example model with a small number of equations.
It consists of one denitrification and 2 nitrification tanks and a settler.

Change into the directory ../ASM1 and translate the model.
Before simulating the model load initial values from the script file small_asm1.mos
that is provided besides the model.
A 14 days dynamic influent data file is provided. So you may simulate up to 14 days.
But start with 1 day as it may take some time for simulation.
After simulation you may have a look at internal concentrations but most interesting
are the relevant concentrations at the effluent of a plant which can be viewed via the
sensors at the effluent of the secondary clarifier.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de
"));
  end SmallPlant;

  class BenchPlant "COST Benchmark WWTP Configuration"
    import WasteWater;

    //Q_air=34574.2654508612 is equal to a Kla of 10 h^-1 from COST benchmark
    //Q_air=12100.99290780142 is equal to a Kla of 3.5 h^-1 from COST benchmark
    extends Modelica.Icons.Example;

    ASM1.EffluentSink Effluent annotation (Placement(transformation(extent={{195,-9},
              {215,11}})));
    ASM1.SludgeSink WasteSludge annotation (Placement(transformation(extent={{211,-44},
              {231,-24}})));
    ASM1.SecClarModTakacs Settler annotation (Placement(transformation(extent={{104,-2},
              {124,18}})));
    ASM1.nitri tank4(V=1333) annotation (Placement(transformation(extent={{-6,-5},
              {14,15}})));
    ASM1.nitri tank3(V=1333) annotation (Placement(transformation(extent={{-30,-5},
              {-10,15}})));
    ASM1.deni tank2 annotation (Placement(transformation(extent={{-77,24},{-57,
              44}})));
    ASM1.deni tank1 annotation (Placement(transformation(extent={{-102,23},{-82,
              43}})));
    ASM1.WWSource WWSource annotation (Placement(transformation(extent={{-114,94},
              {-94,114}})));
    ASM1.sensor_NO sensor_NO annotation (Placement(transformation(extent={{-72,99},
              {-52,119}})));
    ASM1.sensor_O2 sensor_O2 annotation (Placement(transformation(extent={{38,31},
              {55,48}})));
    Modelica.Blocks.Sources.Constant Constant1(k=55338)
                                               annotation (Placement(
          transformation(extent={{-106,-75},{-86,-55}})));
    ASM1.pump WastePump(Q_max=385) annotation (Placement(transformation(extent={{181,-61},
              {201,-41}})));
    Modelica.Blocks.Sources.Constant Constant2(k=18446)
                                               annotation (Placement(
          transformation(extent={{24,-78},{44,-58}})));
    sensor_NH sensor_NH1 annotation (Placement(transformation(extent={{122,31},
              {138,47}})));
    WasteWater.ASM1.sensor_NO sensor_NO1 annotation (Placement(transformation(
            extent={{139,27},{155,43}})));
    WasteWater.ASM1.sensor_TKN sensor_TKN1 annotation (Placement(transformation(
            extent={{152,18},{168,34}})));
    WasteWater.ASM1.sensor_COD sensor_COD1 annotation (Placement(transformation(
            extent={{167,12},{183,28}})));
    WasteWater.ASM1.pump pump annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-90,-14})));
    WasteWater.ASM1.pump pump1 annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={27,-37})));
    WasteWater.ASM1.divider2 divider2_1
      annotation (Placement(transformation(extent={{61,-4},{81,16}})));
    WasteWater.ASM1.mixer3 mixer3_1
      annotation (Placement(transformation(extent={{-135,23},{-115,43}})));
    WasteWater.ASM1.nitri_tank5 tank5
      annotation (Placement(transformation(extent={{19,-5},{39,15}})));
    Modelica.Blocks.Sources.CombiTimeTable CombiTableTime(
      fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/Inf_dry.txt"),
      table=[0,0; 1,1],
      columns=integer(({2,3,4,5,6,7,8,9,10,11,12,13,14,15})),
      tableName="Inf_dry",
      tableOnFile=("Inf_dry") <> "NoName") annotation (Placement(transformation(
            extent={{-143,93},{-123,113}})));
    WasteWater.ASM1.Aeration_energy aeration_energy
      annotation (Placement(transformation(extent={{39,51},{59,71}})));
    WasteWater.ASM1.sensor_Q_EFF sensor_Q_EFF
      annotation (Placement(transformation(extent={{152,-18},{172,2}})));
    WasteWater.ASM1.sensor_SP sensor_SP
      annotation (Placement(transformation(extent={{154,-67},{174,-47}})));
    WasteWater.ASM1.sensor_Q_INF sensor_Q_INF
      annotation (Placement(transformation(extent={{-141,56},{-121,76}})));
    WasteWater.ASM1.sensor_Q sensor_Q
      annotation (Placement(transformation(extent={{77,-24},{57,-4}})));
    WasteWater.ASM1.sensor_Q sensor_Q1
      annotation (Placement(transformation(extent={{105,-40},{85,-20}})));
    WasteWater.ASM1.PE pE
      annotation (Placement(transformation(extent={{129,-91},{149,-71}})));
    WasteWater.ASM1.ME mE
      annotation (Placement(transformation(extent={{40,87},{60,107}})));
    WasteWater.ASM1.sensor_Q sensor_Q2
      annotation (Placement(transformation(extent={{125,-65},{145,-45}})));
    WasteWater.ASM1.sensor_COD sensor_COD
      annotation (Placement(transformation(extent={{-97,59},{-77,79}})));
    WasteWater.ASM1.sensor_COD sensor_COD2
      annotation (Placement(transformation(extent={{-56,53},{-36,73}})));
    WasteWater.ASM1.sensor_COD sensor_COD3
      annotation (Placement(transformation(extent={{-22,24},{-2,44}})));
    WasteWater.ASM1.sensor_COD sensor_COD4
      annotation (Placement(transformation(extent={{0,23},{20,43}})));
    WasteWater.ASM1.sensor_COD sensor_COD5
      annotation (Placement(transformation(extent={{54,30},{74,50}})));
    WasteWater.ASM1.EC eC
      annotation (Placement(transformation(extent={{75,74},{95,94}})));
    WasteWater.ASM1.OCI oCI
      annotation (Placement(transformation(extent={{147,76},{167,96}})));
  equation
    connect(tank3.Out, tank4.In) annotation (Line(points={{-10,5},{-10,5},{-6,5}}));
    connect(tank3.In, tank2.Out) annotation (Line(points={{-30,5},{-44,5},{-44,
            22},{-44,34},{-57,34}}));
    connect(tank1.Out, tank2.In) annotation (Line(points={{-82,33},{-80,33},{
            -80,34},{-77,34}}));
    connect(sensor_NO.In, tank2.MeasurePort) annotation (Line(points={{-62,99},
            {-62,38.5},{-61.5,38.5}}));
    connect(WastePump.Out, WasteSludge.In) annotation (Line(points={{201,-48.2},
            {198,-48.2},{198,-50},{200,-50},{200,-35.2},{211,-35.2}}));
    connect(WastePump.u, Constant2.y)
      annotation (Line(points={{182.1,-48.5},{45,-48.5},{45,-68}},   color={0,0,
            255}));
    connect(sensor_NH1.In, Settler.Effluent) annotation (Line(points={{130,31},
            {130,13.7},{124.2,13.7}}));
    connect(sensor_NO1.In, Settler.Effluent) annotation (Line(points={{147,27},
            {147,13.7},{124.2,13.7}}));
    connect(sensor_TKN1.In, Settler.Effluent) annotation (Line(points={{160,18},
            {160,13.7},{124.2,13.7}}));
    connect(sensor_COD1.In, Settler.Effluent) annotation (Line(points={{175,12},
            {175,13.7},{124.2,13.7}}));

    connect(divider2_1.Out1, Settler.Feed) annotation (Line(
        points={{81,8.6},{101,8.6},{101,9.4},{104,9.4}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pump1.u, Constant2.y) annotation (Line(
        points={{35.9,-39.5},{45,-39.5},{45,-68}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(pump1.Out, mixer3_1.In2) annotation (Line(
        points={{17,-39.8},{-145.5,-39.8},{-145.5,32.5},{-135,32.5}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(mixer3_1.In3, pump.Out) annotation (Line(
        points={{-135,28.5},{-138,28.5},{-138,-16.8},{-100,-16.8}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(mixer3_1.Out, tank1.In) annotation (Line(
        points={{-115,32.6},{-106.5,32.6},{-106.5,33},{-102,33}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(Constant1.y, pump.u) annotation (Line(
        points={{-85,-65},{-65,-65},{-65,-16.5},{-81.1,-16.5}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tank5.In, tank4.Out) annotation (Line(
        points={{19,5},{14,5}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(tank5.Out, divider2_1.In) annotation (Line(
        points={{39,5},{49,5},{49,6.3},{61,6.3}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sensor_O2.In, tank5.MeasurePort) annotation (Line(
        points={{46.5,31},{35,31},{35,9.5},{34.5,9.5}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(CombiTableTime.y, WWSource.data) annotation (Line(
        points={{-122,103},{-114,103},{-114,104},{-113,104}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tank5.Kla, aeration_energy.Kla5) annotation (Line(
        points={{29,8},{29,57},{39,57}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tank4.Kla, aeration_energy.Kla4) annotation (Line(
        points={{4.2,8},{3,8},{3,61},{39,61}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tank3.Kla, aeration_energy.Kla3) annotation (Line(
        points={{-19.8,8},{-20,8},{-20,65},{39,65}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_Q_EFF.In, Settler.Effluent) annotation (Line(
        points={{152,-8},{127,-8},{127,13.7},{124.2,13.7}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_Q_EFF.Out, Effluent.In) annotation (Line(
        points={{172,-8},{185,-8},{185,3},{195,3}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sensor_SP.Out, WastePump.In) annotation (Line(
        points={{174,-57},{174,-54.3},{181,-54.3}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sensor_Q_INF.In, WWSource.Out) annotation (Line(
        points={{-141,66},{-142,66},{-142,80},{-94.2,80},{-94.2,97}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_Q_INF.Out, mixer3_1.In1) annotation (Line(
        points={{-121,66},{-122,66},{-122,46},{-140,46},{-140,36.5},{-135,36.5}},
        color={0,0,0},
        smooth=Smooth.None));

    connect(sensor_Q.In, divider2_1.Out2) annotation (Line(
        points={{77,-14},{90,-14},{90,4.5},{81,4.5}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_Q.Out, pump.In) annotation (Line(
        points={{57,-14},{-13,-14},{-13,-10.7},{-80,-10.7}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sensor_Q1.Out, pump1.In) annotation (Line(
        points={{85,-30},{62,-30},{62,-33.7},{37,-33.7}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sensor_Q1.In, Settler.Return) annotation (Line(
        points={{105,-30},{105,-16.5},{111,-16.5},{111,-1.6}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_Q.Q, pE.Qa) annotation (Line(
        points={{67,-23.8},{68,-23.8},{68,-77},{129,-77}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_Q1.Q, pE.Qr) annotation (Line(
        points={{95,-39.8},{95,-81},{129,-81}},
        color={0,0,127},
        smooth=Smooth.None));

    connect(mE.Kla3, aeration_energy.Kla3) annotation (Line(
        points={{40,96},{-20,96},{-20,65},{39,65}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(mE.Kla4, aeration_energy.Kla4) annotation (Line(
        points={{40,92},{3,92},{3,61},{39,61}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(mE.Kla5, aeration_energy.Kla5) annotation (Line(
        points={{40,88},{29,88},{29,57},{39,57}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_SP.In, sensor_Q2.Out) annotation (Line(
        points={{154.2,-61.2},{150,-61.2},{150,-55},{145,-55}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_Q2.In, Settler.Waste) annotation (Line(
        points={{125,-55},{121,-55},{121,-1.6},{117,-1.6}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_Q2.Q, pE.Qw) annotation (Line(
        points={{135,-64.8},{118,-64.8},{118,-85},{129,-85}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_SP.u, Settler.Return) annotation (Line(
        points={{153.8,-51.8},{132,-51.8},{132,-52},{111,-52},{111,-1.6}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_COD.In, tank1.MeasurePort) annotation (Line(
        points={{-87,59},{-87,37.5},{-86.5,37.5}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_COD2.In, tank2.MeasurePort) annotation (Line(
        points={{-46,53},{-62,53},{-62,38.5},{-61.5,38.5}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(tank3.MeasurePort, sensor_COD3.In) annotation (Line(
        points={{-14.5,9.5},{-14.5,16.25},{-12,16.25},{-12,24}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sensor_COD4.In, tank4.MeasurePort) annotation (Line(
        points={{10,23},{9,23},{9,9.5},{9.5,9.5}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_COD5.In, sensor_O2.In) annotation (Line(
        points={{64,30},{55,30},{55,31},{46.5,31}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(sensor_COD.COD, eC.COD1) annotation (Line(
        points={{-77.2,69},{-0.6,69},{-0.6,93.4},{75,93.4}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_COD2.COD, eC.COD2) annotation (Line(
        points={{-36.2,63},{19.4,63},{19.4,88.4},{75,88.4}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_COD3.COD, eC.COD3) annotation (Line(
        points={{-2.2,34},{35.9,34},{35.9,84},{75,84}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_COD4.COD, eC.COD4) annotation (Line(
        points={{19.8,33},{46,33},{46,79.6},{75,79.6}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_COD5.COD, eC.COD5) annotation (Line(
        points={{73.8,40},{74,40},{74,75.2},{75,75.2}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(aeration_energy.AE, oCI.AE) annotation (Line(
        points={{59.2,61},{103,61},{103,95.4},{147,95.4}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(pE.PE, oCI.PE) annotation (Line(
        points={{149.2,-81},{108,-81},{108,90.8},{147,90.8}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(mE.ME, oCI.ME) annotation (Line(
        points={{59.8,97},{102,97},{102,86},{147,86}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(eC.EC, oCI.EC) annotation (Line(
        points={{95.2,84},{122,84},{122,81.4},{147,81.4}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sensor_SP.SP, oCI.SP) annotation (Line(
        points={{164,-66.8},{134,-66.8},{134,76.8},{147,76.8}},
        color={0,0,127},
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-150,-110},{220,120}},
          grid={1,1}), graphics),
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
"),   experiment(
        StopTime=14,
        __Dymola_NumberOfIntervals=100,
        __Dymola_Algorithm="lsodar"),
      __Dymola_experimentSetupOutput);
  end BenchPlant;

  model ComplexPlant "Complex ASM1 WWTP"
    import WasteWater;

    extends Modelica.Icons.Example;
    ControlledDivider2 cdivider1 annotation (Placement(transformation(extent={{
              -168,65},{-148,85}})));
    Modelica.Blocks.Sources.Constant Constant2(k=0.8)
      annotation (Placement(transformation(extent={{-178,52},{-168,62}})));
    blower blower1(Q_max=162816)
      "there exist 4 blowers of 4240 Nm3/h each, Q_max adusted according active aerated tanks"
       annotation (Placement(transformation(extent={{145,-16},{165,4}})));
    nitri nitri2(
      V=2772,
      alpha=0.305,
      de=5.24,
      R_air=20) annotation (Placement(transformation(extent={{110,18},{130,38}})));
    deni anaerob(V=1287) annotation (Placement(transformation(extent={{-138,13},
              {-118,33}})));
    deni deni1(V=2772) annotation (Placement(transformation(extent={{-80,14},{
              -60,34}})));
    deni deni3(V=2772) annotation (Placement(transformation(extent={{80,18},{
              100,38}})));
    deni deni2(V=2772) annotation (Placement(transformation(extent={{-20,14},{0,
              34}})));
    nitri nitri3(
      V=5602,
      alpha=0.305,
      de=5.24,
      R_air=21) annotation (Placement(transformation(extent={{144,18},{164,38}})));
    blower blower2(Q_max=81408)
      "there exist 4 blowers of 4240 Nm3/h each, Q_max adjusted according active aerated tanks"
       annotation (Placement(transformation(extent={{111,-13},{131,7}})));
    pump ReturnPump(Q_max=60480) annotation (Placement(transformation(
          origin={-44,-94},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    pump RecyclePump(Q_max=60480) annotation (Placement(transformation(
          origin={10,-37},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    pump WastePump(Q_max=1920) annotation (Placement(transformation(extent={{
              128,-104},{148,-84}})));
    ControlledDivider2 cdivider2 annotation (Placement(transformation(extent={{
              -42,68},{-22,88}})));
    EffluentSink Effluent annotation (Placement(transformation(extent={{170,-72},
              {190,-52}})));
    SludgeSink WasteSludge annotation (Placement(transformation(extent={{170,
              -100},{190,-80}})));
    mixer2 mixer2_1 annotation (Placement(transformation(extent={{-45,16},{-25,
              36}})));
    mixer2 mixer2_2 annotation (Placement(transformation(extent={{50,17},{70,37}})));
    mixer3 mixer3_1 annotation (Placement(transformation(extent={{-107,14},{-87,
              34}})));
    mixer2 mixer2_5 annotation (Placement(transformation(extent={{-140,-15},{
              -120,5}})));
    divider2 divider2_1 annotation (Placement(transformation(
          origin={66,-32},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    ControlledDivider2 cdivider3 annotation (Placement(transformation(
          origin={-122,-87},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    Modelica.Blocks.Sources.Constant Constant4(k=0.5)
      annotation (Placement(transformation(extent={{-148,-80},{-138,-70}})));
    nitri nitri1(
      V=5602,
      alpha=0.305,
      de=5.24,
      R_air=21) annotation (Placement(transformation(extent={{12,14},{32,34}})));
    blower blower3(Q_max=162816)
      "there exist 4 blowers of max 4240 Nm3/h, Q_max adusted according active aerated tanks"
       annotation (Placement(transformation(extent={{13,-23},{33,-3}})));
    Modelica.Blocks.Sources.Constant Constant7(k=0.56)
      annotation (Placement(transformation(extent={{-66,98},{-56,108}})));
    mixer2 mixer2_3 annotation (Placement(transformation(extent={{-168,18},{
              -148,38}})));
    PreClar.preclar3 Preclaryfier(V=1372, n_corr=2.138) annotation (Placement(
          transformation(extent={{-136,68},{-116,88}})));
    FlowSource FlowInput annotation (Placement(transformation(extent={{-176,94},
              {-156,114}})));
    ControlledDivider2 ControlledDivider2_1 annotation (Placement(
          transformation(
          origin={-41,-39},
          extent={{-10,-10},{10,10}},
          rotation=180)));
    Modelica.Blocks.Sources.Constant Constant6 annotation (Placement(
          transformation(extent={{-46,-10},{-36,0}})));
    sensor_NO sensor_NO1 annotation (Placement(transformation(extent={{-14,40},
              {6,60}})));
    sensor_NO sensor_NO3 annotation (Placement(transformation(extent={{90,-56},
              {110,-36}})));
    sensor_NH sensor_NH2 annotation (Placement(transformation(extent={{112,-56},
              {132,-36}})));
    sensor_TSS sensor_TSS1 annotation (Placement(transformation(extent={{-6,-86},
              {14,-66}})));
    Modelica.Blocks.Sources.Constant Temperature(k=11.5)
      annotation (Placement(transformation(extent={{58,44},{78,64}})));
    WasteWater.Misc.RecycleController2 RecycleController1(NO3min=1.5)
      annotation (Placement(transformation(extent={{15,-63},{29,-49}})));
    WasteWater.Misc.ReturnController ReturnController1 annotation (Placement(
          transformation(extent={{-28,-115},{-14,-101}})));
    sensor_Q sensor_Q1 annotation (Placement(transformation(extent={{-194,66},{
              -174,86}})));
    WasteWater.Misc.TwoPoint TwoPoint1(
      on=4.5,
      off=4.0,
      out_on=2.5,
      out_off=1.5) annotation (Placement(transformation(extent={{-8,100},{2,110}})));
    WasteWater.Misc.TwoPoint TwoPoint2(
      on=4.5,
      off=4.0,
      out_on=2.5,
      out_off=1.5) annotation (Placement(transformation(extent={{110,100},{120,
              110}})));
    WasteWater.ASM1.sensor_NH sensor_NH1 annotation (Placement(transformation(
            extent={{142,38},{162,58}})));
    Modelica.Blocks.Math.Feedback Feedback1 annotation (Placement(
          transformation(extent={{12,100},{22,110}})));
    Modelica.Blocks.Math.Feedback Feedback2 annotation (Placement(
          transformation(extent={{127,99},{137,109}})));
    WasteWater.ASM1.sensor_O2 sensor_O2_1 annotation (Placement(transformation(
            extent={{164,38},{184,58}})));
    WasteWater.ASM1.sensor_O2 sensor_O2_2 annotation (Placement(transformation(
            extent={{14,38},{34,58}})));
    WasteWater.ASM1.sensor_COD sensor_COD1 annotation (Placement(transformation(
            extent={{134,-56},{154,-36}})));
    WasteWater.ASM1.sensor_COD sensor_COD2 annotation (Placement(transformation(
            extent={{-118,94},{-98,114}})));
    WasteWater.ASM1.Examples.JenaSecClarModTakacs Settler(hsc=3.46, Asc=3704)
      "The depth is calculated based on V and A of the settler and not the true depth."
       annotation (Placement(transformation(extent={{68,-84},{88,-64}})));
    WasteWater.ASM1.sensor_TKN sensor_TKN1 annotation (Placement(transformation(
            extent={{68,-56},{88,-36}})));
    Modelica.Blocks.Sources.CombiTimeTable CombiTableTime1(
      fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/drysim130303.txt"),
      table=[0,0; 1,1],
      columns=integer(({2})),
      tableName="drysim130303",
      tableOnFile=("drysim130303") <> "NoName")
      annotation (Placement(transformation(extent={{-198,96},{-184,112}})));
    Modelica.Blocks.Sources.CombiTimeTable CombiTableTime2(
      fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/drysim130303.txt"),
      table=[0,0; 1,1],
      columns=integer(({3,7,4,6})),
      tableName="drysim130303",
      tableOnFile=("drysim130303") <> "NoName")
      annotation (Placement(transformation(extent={{-134,96},{-118,112}})));
    Modelica.Blocks.Nonlinear.FixedDelay FixedDelay1(delayTime=1/24/6)
      annotation (Placement(transformation(extent={{-16,-63},{-2,-49}})));
    Modelica.Blocks.Math.Feedback Feedback3 annotation (Placement(
          transformation(
          origin={84,-13},
          extent={{5,-5},{-5,5}},
          rotation=180)));
    WasteWater.ASM1.sensor_O2 sensor_O2_3 annotation (Placement(transformation(
            extent={{116,40},{134,58}})));
    Modelica.Blocks.Math.Gain Gain1(k=500) annotation (Placement(transformation(
            extent={{29,95},{49,115}})));
    Modelica.Blocks.Math.Gain Gain2(k=500) annotation (Placement(transformation(
            extent={{144,93},{164,113}})));
    Modelica.Blocks.Math.Gain Gain3(k=500) annotation (Placement(transformation(
            extent={{94,-28},{108,-14}})));
    WasteWater.Misc.TwoPoint TwoPoint3(
      on=4.5,
      off=4.0,
      out_on=2.0,
      out_off=1.0) annotation (Placement(transformation(extent={{63,-19},{75,-7}})));
    Modelica.Blocks.Sources.Step Step1(
      height=0.125,
      offset=-0.5,
      startTime=2.375) annotation (Placement(transformation(extent={{100,-90},{
              110,-80}})));
    WasteWater.ASM1.sensor_TSS sensor_TSS2 annotation (Placement(transformation(
            extent={{166,5},{184,24}})));
  equation
    connect(deni3.Out, nitri2.In) annotation (Line(points={{100,28},{110,28}}));
    connect(nitri2.Out, nitri3.In) annotation (Line(points={{130,28},{144,28}}));
    connect(nitri3.Out, divider2_1.In) annotation (Line(points={{164,28},{190,
            28},{190,-32.3},{76,-32.3}}));
    connect(anaerob.Out, mixer3_1.In2) annotation (Line(points={{-118,23},{
            -110.5,23},{-110.5,23.5},{-107,23.5}}));
    connect(mixer2_5.Out, mixer3_1.In3) annotation (Line(points={{-120,-4.5},{
            -111,-4.5},{-111,19.5},{-107,19.5}}));
    connect(Constant4.y, cdivider3.u)
      annotation (Line(points={{-137.5,-75},{-122,-75},{-122,-81}}, color={0,0,
            255}));
    connect(nitri1.Out, mixer2_2.In2) annotation (Line(points={{32,24},{39,24},
            {39,25.5},{50,25.5}}));
    connect(deni2.Out, nitri1.In) annotation (Line(points={{5.55112e-16,24},{12,
            24}}));
    connect(cdivider2.Out1, mixer2_2.In1) annotation (Line(points={{-22,80.5},{
            40,80.5},{40,29.5},{50,29.5}}));
    connect(cdivider3.Out2, mixer2_5.In1) annotation (Line(points={{-132,-85.5},
            {-172,-85.5},{-172,-2.5},{-140,-2.5}}));
    connect(deni1.Out, mixer2_1.In1) annotation (Line(points={{-60,24},{-56,24},
            {-56,29},{-50.75,29},{-50.75,28.5},{-45,28.5}}));
    connect(Constant7.y, cdivider2.u) annotation (Line(points={{-55.5,103},{
            -50.5,103},{-50.5,54},{-32,54},{-32,72}}, color={0,0,255}));
    connect(cdivider3.In, ReturnPump.Out) annotation (Line(points={{-111.9,
            -87.3},{-86.1,-87.3},{-86.1,-96.8},{-54,-96.8}}));
    connect(cdivider3.Out1, mixer2_3.In2) annotation (Line(points={{-132,-89.5},
            {-184,-89.5},{-184,26.5},{-168,26.5}}));
    connect(cdivider1.Out2, mixer2_3.In1) annotation (Line(points={{-148,73.5},
            {-143,73.5},{-143,40},{-184,40},{-184,30.5},{-168,30.5}}));
    connect(WastePump.Out, WasteSludge.In) annotation (Line(points={{148,-91.2},
            {170,-91.2}}));
    connect(ControlledDivider2_1.Out2, mixer2_1.In2) annotation (Line(points={{-51,
            -37.5},{-52,-37.5},{-52,24.5},{-45,24.5}}));
    connect(ControlledDivider2_1.Out1, mixer2_5.In2) annotation (Line(points={{-51,
            -41.5},{-150,-41.5},{-150,-6.5},{-140,-6.5}}));
    connect(Constant6.y, ControlledDivider2_1.u) annotation (Line(points={{-35.5,
            -5},{-31,-5},{-31,-18},{-41,-18},{-41,-33}}, color={0,0,255}));
    connect(Preclaryfier.In, cdivider1.Out1) annotation (Line(points={{-136,78},
            {-146.5,78},{-146.5,77.5},{-148,77.5}}));
    connect(mixer2_2.Out, deni3.In) annotation (Line(points={{70,27.5},{70,28},
            {80,28}}));
    connect(sensor_NO1.In, deni2.MeasurePort) annotation (Line(points={{-4,40},
            {-4,28.5},{-4.5,28.5}}));
    connect(Temperature.y, anaerob.T) annotation (Line(points={{79,54},{94,54},
            {94,34},{-138,34},{-138,27}}, color={0,0,255}));
    connect(Temperature.y, deni1.T) annotation (Line(points={{79,54},{94,54},{
            94,34},{-80,34},{-80,28}}, color={0,0,255}));
    connect(Temperature.y, deni2.T) annotation (Line(points={{79,54},{94,54},{
            94,34},{-20,34},{-20,28}}, color={0,0,255}));
    connect(Temperature.y, nitri1.T) annotation (Line(points={{79,54},{94,54},{
            94,34},{12,34},{12,28},{12,28}},   color={0,0,255}));
    connect(Temperature.y, deni3.T)
      annotation (Line(points={{79,54},{94,54},{94,34},{80,34},{80,32}}, color=
            {0,0,255}));
    connect(Temperature.y, nitri2.T) annotation (Line(points={{79,54},{94,54},{
            94,34},{110.5,34},{110.5,32},{110,32}}, color={0,0,255}));
    connect(Temperature.y, nitri3.T) annotation (Line(points={{79,54},{94,54},{
            94,34},{144,34},{144,32}}, color={0,0,255}));
    connect(RecycleController1.out, RecyclePump.u) annotation (Line(points={{29.7,
            -56},{34,-56},{34,-39.5},{18.9,-39.5}},      color={0,0,255}));
    connect(FlowInput.Out, sensor_Q1.In) annotation (Line(points={{-156.2,97},{
            -148,97},{-148,90},{-194,90},{-194,76}}));
    connect(sensor_Q1.Q, ReturnController1.in1) annotation (Line(points={{-184,
            66.2},{-184,65.5},{-194,65.5},{-194,-108},{-28.7,-108}},
                                                                   color={0,0,
            255}));
    connect(sensor_NH1.In, nitri3.MeasurePort) annotation (Line(points={{152,38},
            {152,34},{159.5,34},{159.5,32.5}}));
    connect(sensor_O2_2.In, nitri1.MeasurePort) annotation (Line(points={{24,38},
            {24,28.5},{27.5,28.5}}));
    connect(Feedback1.u2, sensor_O2_2.So)
      annotation (Line(points={{17,101},{17,83.5},{33.8,83.5},{33.8,48}},
                                                                      color={0,
            0,255}));
    connect(sensor_O2_1.In, nitri3.MeasurePort) annotation (Line(points={{174,38},
            {174,34},{160.25,34},{160.25,32.5},{159.5,32.5}}));
    connect(Feedback2.u2, sensor_O2_1.So)
      annotation (Line(points={{132,100},{132,81},{183.8,81},{183.8,48}},
                                                                      color={0,
            0,255}));
    connect(sensor_NH1.Snh, TwoPoint1.e) annotation (Line(points={{161.8,48},{
            161.8,72},{-14,72},{-14,105},{-8.5,105}},
                                              color={0,0,255}));
    connect(TwoPoint2.e, sensor_NH1.Snh) annotation (Line(points={{109.5,105},{
            100,105},{100,72},{161.8,72},{161.8,48}},
                                                  color={0,0,255}));
    connect(sensor_COD2.In, Preclaryfier.In) annotation (Line(points={{-108,94},
            {-108,90},{-136,90},{-136,78}}));
    connect(Settler.Effluent, Effluent.In) annotation (Line(points={{88.2,-68.3},
            {108.65,-68.3},{108.65,-68},{128,-68},{128,-60},{170,-60}}));
    connect(sensor_NO3.In, Settler.Effluent) annotation (Line(points={{100,-56},
            {100,-68.3},{88.2,-68.3}}));
    connect(sensor_NH2.In, Settler.Effluent) annotation (Line(points={{122,-56},
            {122,-68.3},{88.2,-68.3}}));
    connect(sensor_COD1.In, Settler.Effluent) annotation (Line(points={{144,-56},
            {144,-60},{128,-60},{128,-68.3},{88.2,-68.3}}));
    connect(WastePump.In, Settler.Waste) annotation (Line(points={{128,-97.3},{
            81,-97.3},{81,-83.6}}));
    connect(Settler.Return, ReturnPump.In) annotation (Line(points={{75,-83.6},
            {75,-90.7},{-34,-90.7}}));
    connect(sensor_TSS1.In, Settler.Return) annotation (Line(points={{4,-86},{4,
            -91},{75,-91},{75,-83.6}}));
    connect(sensor_TKN1.In, Settler.Effluent) annotation (Line(points={{78,-56},
            {88.2,-56},{88.2,-68.3}}));
    connect(CombiTableTime1.y[1], FlowInput.data)
      annotation (Line(points={{-183.3,104},{-175,104}}, color={0,0,255}));
    connect(Preclaryfier.MeasurePort, CombiTableTime2.y) annotation (Line(
          points={{-122.2,87},{-122.2,94},{-116,94},{-116,104},{-117.2,104}},
                                                                          color=
           {0,0,255}));
    connect(FixedDelay1.u, sensor_NO1.Sno) annotation (Line(points={{-17.4,-56},
            {-22,-56},{-22,-8},{10,-8},{10,50},{5.8,50}},
                                                       color={0,0,255}));
    connect(sensor_O2_3.In, nitri2.MeasurePort) annotation (Line(points={{125,40},
            {125,32.5},{125.5,32.5}}));
    connect(sensor_Q1.Out, cdivider1.In) annotation (Line(points={{-174,76},{
            -171,76},{-171,75.3},{-168.1,75.3}}));
    connect(cdivider2.Out2, mixer3_1.In1) annotation (Line(points={{-22,76.5},{
            -18,76.5},{-18,47},{-111,47},{-111,27.5},{-107,27.5}}));
    connect(mixer3_1.Out, deni1.In) annotation (Line(points={{-87,23.6},{-83,
            23.6},{-83,24},{-81,24}}));
    connect(mixer2_1.Out, deni2.In) annotation (Line(points={{-25,26.5},{-22.5,
            26.5},{-22.5,24},{-20,24}}));
    connect(FixedDelay1.y, RecycleController1.in1)
      annotation (Line(points={{-1.3,-56},{14.3,-56}},
                                                     color={0,0,255}));
    connect(Feedback2.u1, TwoPoint2.u)
      annotation (Line(points={{128,104},{124.5,104},{124.5,105},{120,105}},
                                                     color={0,0,255}));
    connect(TwoPoint1.u, Feedback1.u1)
      annotation (Line(points={{2.5,105},{13,105}}, color={0,0,255}));
    connect(mixer2_3.Out, anaerob.In) annotation (Line(points={{-148,28.5},{
            -143,28.5},{-143,23},{-138,23}}));
    connect(Gain1.u, Feedback1.y)
      annotation (Line(points={{27,105},{24,105},{24,105},{21.5,105}},
                                                                     color={0,0,
            255}));
    connect(Gain2.u, Feedback2.y)
      annotation (Line(points={{142,103},{139,103},{139,104},{136.5,104}},
                                                                         color=
            {0,0,255}));
    connect(Feedback3.y, Gain3.u) annotation (Line(points={{88.5,-13},{90.55,
            -13},{90.55,-21},{92.6,-21}}, color={0,0,255}));
    connect(sensor_O2_3.So, Feedback3.u2) annotation (Line(points={{133.82,49},
            {139,49},{139,16},{84,16},{84,-9}},color={0,0,255}));
    connect(TwoPoint3.u, Feedback3.u1)
      annotation (Line(points={{75.6,-13},{80,-13}},
                                                   color={0,0,255}));
    connect(TwoPoint3.e, sensor_NH1.Snh) annotation (Line(points={{62.4,-13},{
            56,-13},{56,6},{104,6},{104,72},{161.8,72},{161.8,48}},
                                                             color={0,0,255}));
    connect(Step1.y, WastePump.u) annotation (Line(points={{110.5,-85},{119,-85},
            {119,-91.5},{129.1,-91.5}},
                                    color={0,0,255}));
    connect(blower1.AirOut, nitri3.AirIn) annotation (Line(points={{154,4},{154,
            5.5},{154,18.2},{154,18.2}}));
    connect(Constant2.y, cdivider1.u)
      annotation (Line(points={{-167.5,57},{-158,57},{-158,69}}, color={0,0,255}));
    connect(Gain1.y, blower3.u) annotation (Line(points={{50,105},{55,105},{55,
            83},{49,83},{49,-16},{32.8,-16}},
                                            color={0,0,255}));
    connect(Gain3.y, blower2.u) annotation (Line(points={{108.7,-21},{139,-21},
            {139,-6},{130.8,-6}}, color={0,0,255}));
    connect(blower1.u, Gain2.y)
      annotation (Line(points={{164.8,-9},{195,-9},{195,103},{165,103}}, color=
            {0,0,255}));
    connect(Settler.Feed, divider2_1.Out1) annotation (Line(points={{68,-72.6},
            {48,-72.6},{48,-34.6},{56,-34.6}}));
    connect(ReturnPump.u, ReturnController1.out) annotation (Line(points={{-35.1,
            -96.5},{-5,-96.5},{-5,-108},{-13.3,-108}}, color={0,0,255}));
    connect(RecyclePump.Out, ControlledDivider2_1.In) annotation (Line(points={{
            -8.88178e-016,-39.8},{-8,-39.8},{-8,-40},{-15.9,-40},{-15.9,-40},{
            -30,-40}}));
    connect(RecyclePump.In, divider2_1.Out2) annotation (Line(points={{20,-33.7},
            {38,-33.7},{38,-30.5},{56,-30.5}}));
    connect(blower3.AirOut, nitri1.AirIn) annotation (Line(points={{22,-3},{22,
            14.2}}));
    connect(cdivider2.In, Preclaryfier.Out) annotation (Line(points={{-42.1,
            78.3},{-80,78.3},{-80,78},{-115,78}}));
    connect(blower2.AirOut, nitri2.AirIn) annotation (Line(points={{120,7},{120,
            18.2}}));
    connect(sensor_TSS2.In, nitri3.Out) annotation (Line(points={{175,5},{175,
            28},{164,28}}));

    annotation (
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-200,-120},{200,120}},
          grid={1,1}), graphics),
      Documentation(info="This ASM1 example plant configuration is from a real municipal wastewater treatment plant
with a size of 145.000 p.e. It is a cascade-type continuous flow plant for a mean dry
weather inflow of 28.500 m3/d. It consists of a preclarifier, an anaerobic tank,
3 denitrification and 3 nitrification tanks and a secondary settler.
This model is an example for the Wastewater library and is not adapted with its parameters
to the reality, therefore simulation results do not reflect the real plant behaviour.

Change into the directory ../ASM1 and translate the model.
Before simulating the model load initial values from the script file complex_asm1.mos
that is provided besides the model.
A 14 days dynamic influent data file is provided. So you may simulate up to 14 days.
But start with 1 day as it may take some time for simulation.
After simulation you may have a look at internal concentrations but most interesting
are the relevant concentrations at the effluent of a plant which can be viewed via the
sensors at the effluent of the secondary clarifier.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de
"));
  end ComplexPlant;

  model JenaSecClarModTakacs
    "Secondary Clarifier Model based on Takacs prepared for a special plant"

    extends WasteWater.Icons.SecClar;
    extends ASM1.SecClar.Takacs.Interfaces.ratios;
    package SCP = ASM1.SecClar.Takacs;
    import SI = Modelica.SIunits;
    package WI = ASM1.Interfaces;
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
