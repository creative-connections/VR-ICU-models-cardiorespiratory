within ;
package modelECMORespiratoryVR

  model VolumeController
    import Modelica.Units.SI.*;
    parameter Volume peakvolume = 0.0005; //peak volume
    parameter Real Iratio = 1 "I in I:E ratio";
    parameter Real Eratio = 2 "E in I:E ratio";
    parameter Real pause = 0 "fraction of RC pause"; //0.01 - 0.75
    parameter Physiolibrary.Types.VolumeFlowRate maxflowrate = 1; //technical limit of ventilator
    parameter Boolean connected = true;

    Physiolibrary.Types.RealIO.VolumeInput volume annotation (Placement(
          transformation(extent={{80,-18},{120,22}}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={76,0}))); //measuredvolume
    Physiolibrary.Types.RealIO.VolumeFlowRateOutput volumeflowrate annotation (
        Placement(transformation(extent={{-60,-76},{-40,-56}}),
          iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,-86})));
    Physiolibrary.Types.RealIO.FrequencyInput frequency annotation (Placement(
          transformation(extent={{-72,54},{-32,94}}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-62,60})));
    Modelica.Blocks.Interfaces.RealOutput outflowopen annotation (Placement(
          transformation(extent={{-88,-62},{-68,-42}}), iconTransformation(
          extent={{-17,-17},{17,17}},
          rotation=180,
          origin={-97,-1})));

      discrete Time T0 "beginning of respiratory cycle";
      Boolean b(start=false);

      discrete Time TE "duration of expiration"; //HP
      discrete Time TP "duration of pause";
      discrete Time TI "duration of inspiration";
      discrete Time RC "duration of respiratory cycle";
      discrete Physiolibrary.Types.VolumeFlowRate inspiratoryflow;
      discrete Physiolibrary.Types.VolumeFlowRate meanventilation2;

      //parameter Frequency HR=1.2;

      Time tr "relative time in respiratory cycle"; //tc

      //parameter Time TD1=0.07 "relative time of start of systole";
      //discrete Time TD2 "relative time of end of systole";
      //parameter MassFlowRate QP=0.424 "peak mass flowrate";
    Physiolibrary.Types.RealIO.VolumeFlowRateOutput meanventilation annotation (
        Placement(transformation(extent={{-60,-76},{-40,-56}}),
          iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={-100,36})));
  initial equation
        T0 = 0 "set beginning of respiratory cycle";
        RC = 1/frequency "update length of respiratory cycle";
        TP = RC * max(min(0.5, pause),0) "compute pause time";
        TI = (RC - TP)* Iratio/(Iratio+Eratio) "compute inspiration duration";
        TE = (RC - TP)* Eratio/(Iratio+Eratio) "compute start of exp phase";
        inspiratoryflow = peakvolume / TI "compute flowrate needed to achieve volume";
        //TD2 = TD1 + (2/5)/HR "compute end time of systole";
         meanventilation2 = peakvolume * frequency;
  equation
      b = time - pre(T0) >= pre(RC) "true if new respiratory cycle begins";
      when {b} then
        //update everything in new respiration cycle
        T0 = time "set beginning of cardiac cycle";
        //the rest same as initial equation
        RC = 1/frequency "update length of respiratory cycle";
        TP = RC * max(min(0.5, pause),0) "compute pause time";
        TI = (RC - TP)* Iratio/(Iratio+Eratio) "compute inspiration duration";
        TE = (RC - TP)* Eratio/(Iratio+Eratio) "compute start of exp phase";
        inspiratoryflow = peakvolume / TI "compute flowrate needed to achieve volume";
        meanventilation2 = peakvolume * frequency;
      end when;
      tr = time - T0 "relative time in respiratory cycle";
      if connected then
      volumeflowrate =
        if tr < TI then inspiratoryflow
        else 0;
      outflowopen =
        if tr < (TI+TP) then 0
        else 1;
        meanventilation = meanventilation2;
      else
        volumeflowrate = 0;
        outflowopen = 1;
        meanventilation = 0;
     end if;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{-80,60},{80,-66}},
            lineColor={28,108,200},
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid), Text(
            extent={{-218,-64},{416,-126}},
            textColor={28,108,200},
            textString="%name")}),            Diagram(coordinateSystem(
            preserveAspectRatio=false)),
      experiment(
        StopTime=90,
        Tolerance=1e-07,
        __Dymola_Algorithm="Dassl"));
  end VolumeController;

  package BloodGasesTransport
    "Transport of O2 and CO2 through respiration and circulation in human body"

    model TissueUnit
      extends Physiolibrary.Icons.SystemicCirculation;
      parameter Physiolibrary.Types.MolarFlowRate O2_consumption=
          1.666666666666667e-05*(2*7.71)
        "Tissue consumption of O2 by metabolism";
      parameter Physiolibrary.Types.MolarFlowRate CO2_production=
          1.666666666666667e-05*(2*6.17)
        "Tissue production of CO2 by metabolism";
      parameter Physiolibrary.Types.HydraulicConductance Conductance=
          1.250102626409427e-07*(1/20) "Tissue blood flow conductance";
      parameter Physiolibrary.Types.Fraction ArteriesViensResistanceRatio=7/8
        "Ratio between arteries and veins resistance";
      parameter Physiolibrary.Types.Volume bloodVolume_start=0.0003;
      parameter Physiolibrary.Types.Volume bloodV0=0.0002;
      parameter Physiolibrary.Types.MassFraction BloodComposition[Blood.nS]=
          Blood.VenousDefault "Initial composition of blood in tissue";
      parameter Physiolibrary.Types.HydraulicCompliance Compliance=
          3.0002463033826e-08 "Compliance of tissue blood vessels";
      Physiolibrary.Fluid.Components.Resistor systemicArteriesResistance(
          redeclare package Medium = Blood, Resistance=1/Conductance*
            ArteriesViensResistanceRatio) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={28,40})));
      Physiolibrary.Fluid.Components.Resistor systemicVeinsResistance(
          redeclare package Medium = Blood, Resistance=1/Conductance*(1 -
            ArteriesViensResistanceRatio)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-64,0})));
      Physiolibrary.Fluid.Components.ElasticVessel systemicCapillaries(
        redeclare package Medium = Blood,
        massFractions_start=BloodComposition,
        useSubstances=true,
        volume_start(displayUnit="l") = bloodVolume_start,
        Compliance(displayUnit="ml/mmHg") = Compliance,
        ZeroPressureVolume(displayUnit="l") = bloodV0,
        nPorts=2) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={-2,-14})));
      Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit = "mmol/min") = O2_consumption) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-70, -38})));
      Chemical.Sources.SubstanceInflowT CO2_left(SubstanceFlow(displayUnit = "mmol/min") = CO2_production, redeclare
          package stateOfMatter =
          Chemical.Interfaces.IdealGas,                                                                                                                                    substanceData = Chemical.Substances.CarbonDioxide_gas()) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {70, -38})));
      Physiolibrary.Fluid.Sensors.BloodGasesMeasurement tissue annotation (
          Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=180,
            origin={38,-8})));
      replaceable package Blood = Physiolibrary.Media.Blood     annotation (
        choicesAllMatching = True);
      Physiolibrary.Fluid.Interfaces.FluidPort_a q_in(redeclare package Medium =
            Blood) annotation (Placement(transformation(rotation=0, extent={{83,
                -5},{96,6}}), iconTransformation(extent={{95,-7},{108,4}})));
      Physiolibrary.Fluid.Interfaces.FluidPort_b q_out(redeclare package Medium =
            Blood) annotation (Placement(transformation(rotation=0, extent={{-97,
                -5},{-84,6}}), iconTransformation(extent={{-105,-5},{-92,6}})));
      Physiolibrary.Types.RealIO.PressureOutput pressure annotation (Placement(
            transformation(extent={{56,-36},{76,-16}}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={80,-64})));
      Physiolibrary.Types.RealIO.PressureOutput pO2 annotation (Placement(
            transformation(extent={{56,4},{76,24}}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-40,-64})));
      Physiolibrary.Types.RealIO.PressureOutput pCO2 annotation (Placement(
            transformation(extent={{56,-8},{76,12}}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-64})));
      Physiolibrary.Types.RealIO.pHOutput pH annotation (Placement(
            transformation(extent={{56,-22},{76,-2}}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={40,-64})));
      Physiolibrary.Types.RealIO.FractionOutput sO2
        "Oxygen saturation (amount of oxygen per amount of hemoglobin units)"
        annotation (Placement(transformation(extent={{56,16},{76,36}}),
            iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-80,-64})));
    equation
      connect(q_in, systemicArteriesResistance.q_in) annotation (
        Line(points = {{89.5, 0.5}, {80, 0.5}, {80, 40}, {38, 40}}, color = {127, 0, 0}));
      connect(q_out, systemicVeinsResistance.q_out) annotation (
        Line(points = {{-90.5, 0.5}, {-74, 7.77156e-16}}, color = {127, 0, 0}));
      connect(tissue.pressure, pressure) annotation (
        Line(points = {{49, -16}, {49, -26}, {66, -26}}, color = {0, 0, 127}));
      connect(tissue.pO2, pO2) annotation (
        Line(points = {{49, -4}, {52, -4}, {52, 14}, {66, 14}}, color = {0, 0, 127}));
      connect(tissue.pCO2, pCO2) annotation (
        Line(points = {{49, -8}, {54, -8}, {54, 2}, {66, 2}}, color = {0, 0, 127}));
      connect(tissue.pH, pH) annotation (
        Line(points = {{49, -12}, {52, -12}, {52, -12}, {66, -12}}, color = {0, 0, 127}));
      connect(systemicArteriesResistance.q_out, systemicCapillaries.q_in[1]) annotation (
        Line(points={{18,40},{-2,40},{-2,-12},{-1.9,-12},{-1.9,-14.65}},           color = {127, 0, 0}, thickness = 0.5));
      connect(systemicCapillaries.q_in[2], tissue.a_port) annotation (
        Line(points={{-1.9,-13.35},{-2,-13.35},{-2,-28},{32,-28},{32,-18.2}},          color = {127, 0, 0}, thickness = 0.5));
      connect(tissue.b_port, systemicVeinsResistance.q_in) annotation (
        Line(points = {{44.2, -18.2}, {44.2, -32}, {-48, -32}, {-48, -1.66533e-15}, {-54, -1.66533e-15}}, color = {127, 0, 0}, thickness = 0.5));
      connect(tissue.sO2, sO2) annotation (
        Line(points = {{49, -3.55271e-15}, {50, -3.55271e-15}, {50, 26}, {66, 26}}, color = {0, 0, 127}));
      connect(O2_left.port_a, systemicCapillaries.substances.O2) annotation (
        Line(points = {{-60, -38}, {14, -38}, {14, -14}, {8, -14}}, color = {158, 66, 200}));
      connect(systemicCapillaries.substances.CO2, CO2_left.port_b) annotation (
        Line(points = {{8, -14}, {14, -14}, {14, -38}, {60, -38}}, color = {158, 66, 200}));
      connect(tissue.substances, systemicCapillaries.substances)
        annotation (Line(
          points={{28,-8},{16,-8},{16,-14},{8,-14}},
          color={158,66,200},
          thickness=0.5));
      annotation (
        Diagram(coordinateSystem(extent = {{-90, -50}, {90, 50}})),
        Icon(coordinateSystem(extent = {{-90, -50}, {90, 50}})));
    end TissueUnit;

    model RespiratoryCenter
      extends Physiolibrary.Icons.RespiratoryCenter;
      Modelica.Blocks.Math.Add add annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 90, origin = {-59, 15})));
      Physiolibrary.Types.Constants.PressureConst pressure(k(displayUnit = "kPa") = -pCO2_ZeroVentilation) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin = {-56, -10})));
      Physiolibrary.Types.Constants.HydraulicConductanceConst hydraulicConductance2(k(displayUnit = "ml/(kPa.min)") = 2.5e-07) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin = {-68, -48})));
      Modelica.Blocks.Math.Product product3 annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-52, 40})));
      Modelica.Blocks.Math.Max totalVentilation annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 90, origin = {-67, 69})));
      Modelica.Blocks.Sources.Constant const(k = 0) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin = {-72, 52})));
      Modelica.Blocks.Math.Add add1(k2 = +1) annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin = {52, -2})));
      Physiolibrary.Types.Constants.VolumeConst deadspace_c(k(displayUnit = "l") = 0.00035) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin = {70, -38})));
      Modelica.Blocks.Math.Division respirationRate annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin = {50, 42})));
      Modelica.Blocks.Math.Product RR_multiply annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 180, origin = {39, 63})));
      Physiolibrary.Types.Constants.VolumeConst base_deadspace_volume(k = DV) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 180, origin = {70, 64})));
      Physiolibrary.Types.Constants.VolumeConst volume(k(displayUnit = "l") = 0.0023) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin = {68, -2})));
      Modelica.Blocks.Math.Min tidalVolume annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 90, origin = {55, 21})));
      Modelica.Blocks.Math.Division slope annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin = {-30, 0})));
      Modelica.Blocks.Math.Max max1 annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 90, origin = {-27, -29})));
      Modelica.Blocks.Math.Add add2(k2 = -1) annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin = {-40, -58})));
      Modelica.Blocks.Sources.Constant const2(k = 5000 - 4300) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 90, origin = {-23, -63})));
      Modelica.Blocks.Math.Product product4 annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin = {36, -38})));
      Modelica.Blocks.Sources.Constant const3(k = 4300) annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin = {-36, -88})));
      Modelica.Blocks.Math.Gain W(k = 90 * (101325 / 760) - 4.3 * 1000) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 90, origin = {-59, -31})));
      Physiolibrary.Types.Constants.FrequencyConst m(k = 0.505) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 180, origin = {70, -70})));
      parameter Physiolibrary.Types.Pressure pCO2_ZeroVentilation(displayUnit = "kPa") = 4800 "Long term adaptation for pCO2 during acidosis/alcalosis";
      parameter Physiolibrary.Types.Volume DV = 0.00015 "Deadspace volume";
      Modelica.Blocks.Interfaces.RealOutput deadSpaceVentilation annotation (
        Placement(transformation(rotation = 0, extent = {{98, 21.5}, {118, 46.5}}), iconTransformation(extent = {{-10, -12.5}, {10, 12.5}}, rotation = 270, origin = {60, -108})));
      Modelica.Blocks.Interfaces.RealInput pCO2 "Arterial pCO2" annotation (
        Placement(transformation(rotation = 0, extent = {{-96, -20.5}, {-76, 4.5}}), iconTransformation(extent = {{-10, -12.5}, {10, 12.5}}, rotation = 0, origin = {-104, -10})));
      Modelica.Blocks.Interfaces.RealInput pO2 "Arterial pO2" annotation (
        Placement(transformation(rotation = 0, extent = {{-82, -96.5}, {-62, -71.5}}), iconTransformation(extent = {{-10, -12.5}, {10, 12.5}}, rotation = 0, origin = {-106, -90})));
      Modelica.Blocks.Interfaces.RealOutput ventilation annotation (
        Placement(transformation(rotation = 0, extent={{2,74.5},{22,99.5}}),         iconTransformation(extent = {{-10, -12.5}, {10, 12.5}}, rotation = 270, origin = {4, -107})));
      Modelica.Blocks.Logical.Switch VentilationSwitch
        annotation (Placement(transformation(extent={{-22,74},{-2,94}})));
      Modelica.Blocks.Sources.BooleanConstant physiologyVentilation
        "if true then computed ventilation rate, false - constant parameter may be driven from outside"
        annotation (Placement(transformation(extent={{-46,78},{-36,88}})));
      Physiolibrary.Types.Constants.VolumeFlowRateConst ArtificialVentilation(k
          =8.3333333333333e-05)
        annotation (Placement(transformation(extent={{-36,50},{-28,58}})));
      Modelica.Blocks.Continuous.Filter filter(
        analogFilter=Modelica.Blocks.Types.AnalogFilter.Butterworth,
        filterType=Modelica.Blocks.Types.FilterType.LowPass,
        order=4,
        f_cut=0.5)
        annotation (Placement(transformation(extent={{-20,38},{-8,50}})));
    equation
      connect(const.y, totalVentilation.u1) annotation (
        Line(points = {{-72, 56.4}, {-71.2, 60.6}}, color = {0, 0, 127}));
      connect(pressure.y, add.u2) annotation (
        Line(points = {{-56, -5}, {-46, -5}, {-46, 6.6}, {-54.8, 6.6}}, color = {0, 0, 127}));
      connect(add.y, product3.u1) annotation (
        Line(points = {{-59, 22.7}, {-72, 22.7}, {-72, 28}, {-58, 28}}, color = {0, 0, 127}));
      connect(totalVentilation.u2, product3.y) annotation (
        Line(points = {{-62.8, 60.6}, {-52, 60.6}, {-52, 51}}, color = {0, 0, 127}));
      connect(deadspace_c.y, add1.u2) annotation (
        Line(points = {{70, -33}, {70, -26}, {55.6, -26}, {55.6, -9.2}}, color = {0, 0, 127}));
      connect(base_deadspace_volume.y, RR_multiply.u2) annotation (
        Line(points = {{65, 64}, {65, 67.2}, {47.4, 67.2}}, color = {0, 0, 127}));
      connect(respirationRate.y, RR_multiply.u1) annotation (
        Line(points = {{50, 48.6}, {62, 48.6}, {62, 56}, {54, 56}, {54, 58.8}, {47.4, 58.8}}, color = {0, 0, 127}));
      connect(tidalVolume.u1, add1.y) annotation (
        Line(points = {{52, 15}, {52, 4.6}}, color = {0, 0, 127}));
      connect(tidalVolume.y, respirationRate.u2) annotation (
        Line(points = {{55, 26.5}, {53.6, 26.5}, {53.6, 34.8}}, color = {0, 0, 127}));
      connect(volume.y, tidalVolume.u2) annotation (
        Line(points = {{68, 3}, {68, 15}, {58, 15}}, color = {0, 0, 127}));
      connect(const2.y, max1.u2) annotation (
        Line(points = {{-23, -57.5}, {-26, -57.5}, {-26, -42}, {-24, -42}, {-24, -35}}, color = {0, 0, 127}));
      connect(add2.y, max1.u1) annotation (
        Line(points = {{-40, -51.4}, {-30, -51.4}, {-30, -35}}, color = {0, 0, 127}));
      connect(max1.y, slope.u2) annotation (
        Line(points = {{-27, -23.5}, {-26, -23.5}, {-26, -7.2}, {-26.4, -7.2}}, color = {0, 0, 127}));
      connect(const3.y, add2.u2) annotation (
        Line(points = {{-36, -81.4}, {-32, -81.4}, {-32, -78}, {-36, -78}, {-36, -65.2}, {-36.4, -65.2}}, color = {0, 0, 127}));
      connect(product4.y, add1.u1) annotation (
        Line(points = {{36, -31.4}, {36, -9.2}, {48.4, -9.2}}, color = {0, 0, 127}));
      connect(slope.y, product3.u2) annotation (
        Line(points = {{-30, 6.6}, {-30, 18}, {-46, 18}, {-46, 28}}, color = {0, 0, 127}));
      connect(hydraulicConductance2.y, W.u) annotation (
        Line(points = {{-68, -43}, {-68, -37}, {-59, -37}}, color = {0, 0, 127}));
      connect(W.y, slope.u1) annotation (
        Line(points = {{-59, -25.5}, {-59, -18}, {-33.6, -18}, {-33.6, -7.2}}, color = {0, 0, 127}));
      connect(m.y, product4.u2) annotation (
        Line(points = {{65, -70}, {48, -70}, {48, -64}, {39.6, -64}, {39.6, -45.2}}, color = {0, 0, 127}));
      connect(deadSpaceVentilation, RR_multiply.y) annotation (
        Line(points = {{108, 34}, {90, 34}, {90, 90}, {26, 90}, {26, 63}, {31.3, 63}}, color = {0, 0, 127}));
      connect(pCO2, add.u1) annotation (
        Line(points = {{-86, -8}, {-63.2, -8}, {-63.2, 6.6}}, color = {0, 0, 127}));
      connect(pO2, add2.u1) annotation (
        Line(points = {{-72, -84}, {-52, -84}, {-52, -65.2}, {-43.6, -65.2}}, color = {0, 0, 127}));
      connect(physiologyVentilation.y, VentilationSwitch.u2) annotation (Line(
            points={{-35.5,83},{-35.5,84},{-24,84}}, color={255,0,255}));
      connect(totalVentilation.y, VentilationSwitch.u1) annotation (Line(points
            ={{-67,76.7},{-68,76.7},{-68,92},{-24,92}}, color={0,0,127}));
      connect(ventilation, VentilationSwitch.y) annotation (Line(points={{12,87},
              {12,88},{6,88},{6,84},{-1,84}}, color={0,0,127}));
      connect(VentilationSwitch.y, respirationRate.u1) annotation (Line(points=
              {{-1,84},{0,84},{0,34.8},{46.4,34.8}}, color={0,0,127}));
      connect(VentilationSwitch.y, product4.u1) annotation (Line(points={{-1,84},
              {2,84},{2,-56},{32.4,-56},{32.4,-45.2}}, color={0,0,127}));
      connect(ArtificialVentilation.y, filter.u) annotation (Line(points={{-27,
              54},{-24,54},{-24,44},{-21.2,44}}, color={0,0,127}));
      connect(filter.y, VentilationSwitch.u3) annotation (Line(points={{-7.4,44},
              {-8,44},{-8,64},{-24,64},{-24,76}}, color={0,0,127}));
      annotation (
        Diagram(coordinateSystem(extent={{-140,-100},{100,100}}),
                graphics={  Text(extent = {{6, -76}, {34, -84}}, textColor = {28, 108, 200}, textString = "Calculation of slope")}), Icon(
            coordinateSystem(extent={{-140,-100},{100,100}})));
    end RespiratoryCenter;

    model BloodyMaryPPG2
      replaceable package Blood =
        Physiolibrary.Media.Blood                                             annotation (
        choicesAllMatching = True);
      replaceable package Air = Physiolibrary.Media.Air annotation (
        choicesAllMatching = True);
      parameter Physiolibrary.Types.Frequency RR=0.286   "Respiration rate";
      parameter Physiolibrary.Types.Volume TV=0.0005   "Tidal volume";
      parameter Physiolibrary.Types.Volume DV=0.00015   "Dead space volume";
      parameter Physiolibrary.Types.VolumeFlowRate CO=9.1666666666667e-05   "Cardiac output";
      parameter Real cShuntFrac = 0.02;
      parameter Physiolibrary.Types.HydraulicConductance cShunt=1.250102626409427e-07*(1/3*cShuntFrac);
      //orig 1.250102626409427e-07*(1/3*0.02);
      parameter Physiolibrary.Types.HydraulicConductance cTotalVentilation=1.019716212977928e-05
          *(1/1.5);
      parameter Physiolibrary.Types.HydraulicConductance cTotalCirculation=1.250102626409427e-07
          *(1/3*(1 - 0.02));
      parameter Physiolibrary.Types.HydraulicCompliance LungsCompliance=6.0004926067653e-07 "Lungs compliance";
      parameter Physiolibrary.Types.Volume ResidualVolume=0.0013 "Lungs residual volume";
      parameter Physiolibrary.Types.Volume TotalCapacity=0.00623 "Lungs total capacity";
      parameter Physiolibrary.Types.Volume BaseTidalVolume=0.0005 "Lungs base tidal volume";
      parameter Physiolibrary.Types.Volume alveolarVolume_start = TV - DV + alveolarV0 "initial volume of air in alveoli";
      parameter Physiolibrary.Types.Volume alveolarV0=0.0013   "volume of air in alveoli, which does not generate air pressure";
      parameter Physiolibrary.Types.Volume lungCapyVolume_start=0.00015   "initial volume of blood in alveolar capillaries";
      parameter Physiolibrary.Types.Volume lungCapyV0=0.0001   "volume of blood in alveolar capillaries, which does not generate blood pressure";
      parameter Physiolibrary.Types.Volume tissueBloodVolume_start=0.0003   "initial volume of blood in tissues";
      parameter Physiolibrary.Types.Volume tissueV0=0.0002   "volume of blood in tissues, which does not generate blood pressure";
      parameter Physiolibrary.Types.HydraulicCompliance CapillariesCompliance=3.0002463033826e-08
                                                                                     "Systemic capillaries compliance";
      parameter Physiolibrary.Types.MolarFlowRate O2_consumption=1.666666666666667e-05
          *(2*7.71) "Tissue consumption of O2 by metabolism";
      parameter Physiolibrary.Types.MolarFlowRate CO2_production=1.666666666666667e-05
          *(2*6.17) "Tissue production of CO2 by metabolism";
      parameter Physiolibrary.Types.HydraulicConductance TotalSystemicConductance=1.250102626409427e-07
          *(1/20) "Total systemic blood circulation conductance";
      parameter Integer NA=1  "Number of pulmonary alveolar units";
      parameter Integer NT=1  "Number of systemic tissue units";
    public
      parameter Physiolibrary.Types.MassFraction ArterialBloodComposition[Blood.nS]
        =Blood.ArterialDefault "Initial composition of arterial blood";
      parameter Physiolibrary.Types.MassFraction VenousBloodComposition[Blood.nS]=
          Blood.VenousDefault "Initial composition of venous blood";
    //  parameter Physiolibrary.Types.Fraction AirO2=0.21 "O2 content in inspired air";
      parameter Physiolibrary.Types.Fraction AirCO2=0.0003
        "CO2 content in inspired air";
      parameter Physiolibrary.Types.Fraction AirH2O=0.06
        "H2O content in inspired air";
      Physiolibrary.Types.Fraction AirN2 "N2 content in inspired air"; //=1 - AirO2.y - AirCO2 - AirH2O
      parameter Physiolibrary.Types.Fraction AirN2_start = 1 - AirO2Fraction.k - AirCO2 - AirH2O;
      Physiolibrary.Types.Pressure Air_pO2=AirO2.y*system.p_ambient
        "O2 content in inspired air";
      Physiolibrary.Types.Pressure Air_pCO2=AirCO2*system.p_ambient
        "CO2 content in inspired air";
      Physiolibrary.Fluid.Components.VolumePump deadSpaceVentilation(redeclare
          package Medium =                                                                      Air, useSolutionFlowInput = true, SolutionFlow = DV * RR,
        density(start=1.1300953420756321, displayUnit="g/cm3"),
        q_in(m_flow(start=5.759563136974551E-05)),
        q_out(p(start=101335.86111397855, displayUnit="bar")))                                                                                            annotation (
        Placement(transformation(extent = {{-14, -52}, {6, -32}})));
      Physiolibrary.Fluid.Sources.PressureSource pressureSource(pressure_start(displayUnit = "Pa"), redeclare
          package Medium =                                                                                                     Air,  massFractions_start = Air.X(AirO2Fraction.k, AirCO2, AirH2O, AirN2_start)) annotation (
        Placement(transformation(extent = {{-96, -52}, {-76, -32}})));
      Physiolibrary.Fluid.Sources.VolumeOutflowSource volumeOutflow(useSolutionFlowInput = true, SolutionFlow = TV * RR, redeclare
          package Medium =                                                                                                                          Air,
        density(start=1.1300953420756321))                                                                                                               annotation (
        Placement(transformation(extent = {{64, -52}, {84, -32}})));
      Physiolibrary.Fluid.Components.VolumePump leftHeartPump(redeclare package
          Medium =                                                                       Blood, useSolutionFlowInput = true, SolutionFlow(displayUnit = "l/min") = 8.3333333333333e-05) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {42, -150})));
      Physiolibrary.Fluid.Components.VolumePump rightHeartPump(redeclare
          package Medium =                                                                Blood, useSolutionFlowInput = true, SolutionFlow(displayUnit = "l/min") = 8.3333333333333e-05) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-50, -150})));
      Physiolibrary.Fluid.Components.ElasticVessel pulmonaryArteries(redeclare
          package Medium =                                                                      Blood,  massFractions_start = VenousBloodComposition, volume_start(displayUnit = "l") = 0.00038, Compliance(displayUnit = "ml/mmHg") = 3.6002955640592e-08, ZeroPressureVolume(displayUnit = "l") = 0.0003, nPorts = 4) annotation (
        Placement(transformation(extent = {{-60, -112}, {-40, -92}})));
      Physiolibrary.Fluid.Components.ElasticVessel pulmonaryVeins(redeclare
          package Medium =                                                                   Blood, massFractions_start = ArterialBloodComposition, volume_start(displayUnit = "l") = 0.0004, Compliance(displayUnit = "ml/mmHg") = 7.5006157584566e-08, ZeroPressureVolume(displayUnit = "l") = 0.0004, nPorts = 4) annotation (
        Placement(transformation(extent = {{32, -112}, {52, -92}})));
      Physiolibrary.Fluid.Components.ElasticVessel systemicArteries(redeclare
          package Medium = Blood,
          massFractions_start = ArterialBloodComposition,
          useSubstances = true,
          volume_start(displayUnit = "l") = 0.00085,
          Compliance(displayUnit = "ml/mmHg") = 2.6627185942521e-08,
          ZeroPressureVolume(displayUnit = "l") = 0.00045,
          nPorts = NT+1,
        chemicalSolution(bloodGases(
            cHCO3(start=24.51879571586536, displayUnit="mmol/l"),
            pCO(start=1.0005449962821231E-05, displayUnit="bar"),
            pCO2(start=5345.09227220704, displayUnit="bar"),
            pO2(start=11608.744631309959, displayUnit="bar"),
            sO2CO(start=0.9627613894820818, displayUnit="%"))))                                                                                                                                                                                                         annotation (
        Placement(transformation(extent = {{46, -206}, {26, -186}})));
      Physiolibrary.Fluid.Components.ElasticVessel systemicVeins(redeclare
          package Medium =
                   Blood,
          massFractions_start = VenousBloodComposition,
          useSubstances = true, volume_start(displayUnit = "l") = 0.00325,
          Compliance(displayUnit = "ml/mmHg") = 6.1880080007267e-07,
          ZeroPressureVolume(displayUnit = "l") = 0.00295,
          nPorts = NT+2,
        chemicalSolution(bloodGases(
            cHCO3(start=26.674314102391666, displayUnit="mmol/l"),
            pCO(start=6.3972838566901375E-06, displayUnit="bar"),
            pCO2(start=6930.575174544387, displayUnit="bar"),
            pO2(start=5006.216473490174, displayUnit="bar"),
            sCO(start=1.808984022893137E-07, displayUnit="%"),
            sO2CO(start=0.649370212847236, displayUnit="%"))))                                                                                                                                                                                                         annotation (
        Placement(transformation(extent = {{-60, -206}, {-40, -186}})));
      Physiolibrary.Fluid.Sensors.PressureMeasure pressureMeasureVeins(redeclare
          package Medium =                                                                        Blood) annotation (
        Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 0, origin = {-80, -204})));
      Physiolibrary.Fluid.Sensors.PressureMeasure pressureMeasurePulmArteries(redeclare
          package Medium =                                                                               Blood) annotation (
        Placement(transformation(extent = {{-60, -108}, {-80, -88}})));
      Physiolibrary.Fluid.Sensors.PressureMeasure pressureMeasurePulmVeins(redeclare
          package Medium =                                                                            Blood) annotation (
        Placement(transformation(extent = {{52, -108}, {72, -88}})));
      Physiolibrary.Types.Constants.VolumeFlowRateConst leftCardiacOutput(k=CO)
        annotation (Placement(transformation(
            extent={{-4,-4},{4,4}},
            rotation=180,
            origin={76,-150})));
      Modelica.Blocks.Math.MultiProduct multiProduct1(nu = 2) annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 0, origin = {-70, -150})));
      Physiolibrary.Types.Constants.HydraulicConductanceConst hydraulicConductance1(k = 1.250102626409427e-07 * (5 / 4)) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 270, origin = {-80, -132})));
      Physiolibrary.Fluid.Components.Conductor pulmonaryShunt(redeclare package
          Medium =                                                                       Blood, Conductance(displayUnit = "l/(cmH2O.s)") = cShunt) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-4, -104})));
    Physiolibrary.Organs.Lungs.Components.RespiratoryUnit respiratoryUnit[NA](
        redeclare package Blood = Blood,
        redeclare package Air = Air,
        AirVolume_initial=ones(NA)*(alveolarVolume_start/NA),
        FunctionalResidualCapacity=ones(NA)*(alveolarV0/NA),
        TotalCompliance=ones(NA)*(LungsCompliance/NA),
        ResidualVolume=ones(NA)*(ResidualVolume/NA),
        TotalCapacity=ones(NA)*(TotalCapacity/NA),
        BaseTidalVolume=ones(NA)*(BaseTidalVolume/NA),
        TotalResistance=ones(NA)*(NA/cTotalVentilation),
        CapillariesVolume_initial=ones(NA)*(lungCapyVolume_start/NA),
        Blood_initial=fill(ArterialBloodComposition, NA),
        CapillariesZeroPressureVolume=ones(NA)*(lungCapyV0/NA),
        CapillariesCompliance=fill(CapillariesCompliance/NA, NA),
        CapillariesConductance=ones(NA)*((cTotalCirculation)/(NA*1/8)),
        ArteriesConductance=ones(NA)*((cTotalCirculation)/(NA*7/8)),
        each nPorts=2,
        pulmCapysVentilated(chemicalSolution(bloodGases(
              cHCO3(start=fill(24.518795715865362, NA), each displayUnit="mmol/l"),
              pCO(start=fill(1.0005449962821248E-05, NA), each displayUnit="bar"),
              pCO2(start=fill(5345.092272207041, NA), each displayUnit="bar"),
              pO2(start=fill(11608.74463130998, NA), each displayUnit="bar"),
              sO2CO(start=fill(0.962761389482082, NA), each displayUnit="%")))))
        annotation (Placement(transformation(rotation=0, extent={{-22,-96},{10,-64}})));

      Physiolibrary.Fluid.Sensors.BloodGasesMeasurement arterial(redeclare
          package Medium =
                   Blood) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=270,
            origin={54,-176})));
      TissueUnit tissueUnit[NT](
        O2_consumption = fill(O2_consumption / NT, NT),
        CO2_production = fill(CO2_production / NT, NT),
        Conductance = fill(TotalSystemicConductance / NT, NT),
        bloodVolume_start = fill(tissueBloodVolume_start / NT, NT),
        bloodV0 = fill(tissueV0 / NT, NT),
        BloodComposition = fill(VenousBloodComposition, NT),
        Compliance=fill(CapillariesCompliance/NT, NT),       redeclare package
          Blood =          Blood,
        systemicCapillaries(chemicalSolution(bloodGases(
              cHCO3(start=fill(26.674314102391666,NT), each displayUnit="mmol/l"),
              pCO(start=fill(6.39728383892192E-06,NT), each displayUnit="bar"),
              pCO2(start=fill(6930.57517454441,NT), each displayUnit="bar"),
              pO2(start=fill(5006.216473490139,NT), each displayUnit="bar"),
              sCO(start=fill(1.8089840228933077E-07,NT), each displayUnit="%"),
              sO2CO(start=fill(0.6493702128472361,NT), each displayUnit="%")))))
         annotation (
        Placement(transformation(extent = {{-14, -202}, {12, -188}})));
        //,
      /*  systemicCapillaries(chemicalSolution(bloodGases(
          cHCO3(start={26.674314102391666}, each displayUnit="mmol/l"),
          pCO(start={6.39728383892192E-06}, each displayUnit="mmHg"),
          pCO2(start={6930.57517454441}, each displayUnit="mmHg"),
          pO2(start={5006.216473490139}, each displayUnit="mmHg"),
          sCO(start={1.8089840228933077E-07}, each displayUnit="%"),
          sO2CO(start={0.6493702128472361}, each displayUnit="%"))))*/
      RespiratoryCenter respiratoryCenter(DV = DV, ArtificialVentilation(k=0.00015),
        physiologyVentilation(k=false),
        VentilationSwitch(y(start=0.00015)))       annotation (
        Placement(transformation(extent={{58,4},{22,34}})));
      inner Modelica.Fluid.System system(T_ambient = 310.15) annotation (
        Placement(transformation(extent={{-94,12},{-74,32}})));
      Physiolibrary.Fluid.Sensors.BloodGasesMeasurement venous(redeclare
          package Medium =
                   Blood) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-66,-176})));
      Physiolibrary.Types.Constants.PressureConst pressure(k(displayUnit="mmHg") = -533.28954966)
        annotation (Placement(transformation(extent={{-42,-58},{-34,-50}})));
      Modelica.Blocks.Interfaces.RealOutput plethy "so2*arterialpressure"
        annotation (Placement(transformation(extent={{94,-204},{114,-184}})));
      Physiolibrary.Types.Constants.PressureConst hemodynamicPressure(k(displayUnit
            ="mmHg") = 120*133.322)
        annotation (Placement(transformation(extent={{68,-212},{76,-204}})));
      Modelica.Blocks.Math.Product PPG
        annotation (Placement(transformation(extent={{78,-198},{90,-186}})));
      Physiolibrary.Types.Constants.FractionConst AirO2Fraction(k=0.21)
        annotation (Placement(transformation(extent={{-62,22},{-54,30}})));
      Modelica.Blocks.Continuous.Filter AirO2(
        analogFilter=Modelica.Blocks.Types.AnalogFilter.Butterworth,
        filterType=Modelica.Blocks.Types.FilterType.LowPass,
        order=4,
        f_cut=0.5)
        annotation (Placement(transformation(extent={{-32,10},{-12,30}})));
    equation
      AirN2=1 - AirO2.y - AirCO2 - AirH2O;
      connect(deadSpaceVentilation.q_out, volumeOutflow.q_in) annotation (
        Line(points = {{6, -42}, {64, -42}}, color = {127, 0, 0}, thickness = 0.5));
      connect(pressureSource.y, deadSpaceVentilation.q_in) annotation (
        Line(points = {{-76, -42}, {-14, -42}}, color = {127, 0, 0}, thickness = 0.5));
      connect(rightHeartPump.q_out, pulmonaryArteries.q_in[1]) annotation (
        Line(points={{-50,-140},{-50,-116},{-50.1,-116},{-50.1,-102.975}},         color = {127, 0, 0}, thickness = 0.5));
      connect(leftHeartPump.q_in, pulmonaryVeins.q_in[1]) annotation (
        Line(points={{42,-140},{42,-100},{41.9,-100},{41.9,-102.975}},         color = {127, 0, 0}, thickness = 0.5));
      connect(pressureMeasureVeins.port, systemicVeins.q_in[1]) annotation (
        Line(points={{-80,-214},{-50,-214},{-50,-196},{-50.1,-196}},                  color = {127, 0, 0}, thickness = 0.5));
      connect(pressureMeasurePulmArteries.port, pulmonaryArteries.q_in[2]) annotation (
        Line(points={{-70,-108},{-50,-108},{-50,-102.325},{-50.1,-102.325}},        color = {127, 0, 0}, thickness = 0.5));
      connect(pulmonaryVeins.q_in[2],pressureMeasurePulmVeins.port)  annotation (
        Line(points={{41.9,-102.325},{42,-102.325},{42,-108},{62,-108}},        color = {127, 0, 0}, thickness = 0.5));
      connect(multiProduct1.y, rightHeartPump.solutionFlow) annotation (
        Line(points = {{-62.98, -150}, {-57, -150}}, color = {0, 0, 127}));
      connect(hydraulicConductance1.y, multiProduct1.u[1]) annotation (
        Line(points={{-80,-137},{-80,-148},{-76,-148},{-76,-151.05}},         color = {0, 0, 127}));
      connect(pressureMeasureVeins.pressure, multiProduct1.u[2]) annotation (
        Line(points={{-86,-208},{-94,-208},{-94,-148.95},{-76,-148.95}},        color = {0, 0, 127}));
      connect(pulmonaryShunt.q_in, pulmonaryArteries.q_in[3]) annotation (
        Line(points={{-14,-104},{-50.1,-104},{-50.1,-101.675}},       color = {127, 0, 0}, thickness = 0.5));
      connect(pulmonaryShunt.q_out, pulmonaryVeins.q_in[3]) annotation (
        Line(points={{6,-104},{41.9,-104},{41.9,-101.675}},       color = {127, 0, 0}, thickness = 0.5));
      for i in 1:NA loop
        connect(respiratoryUnit[i].blood_in, pulmonaryArteries.q_in[4]) annotation (
            Line(
            points={{-22,-87.04},{-30,-87.04},{-30,-100},{-50.1,-100},{-50.1,-101.025}},
            color={127,0,0},
            thickness=0.5));
        connect(respiratoryUnit[i].airways[1], pressureSource.y) annotation (Line(
            points={{-6,-64.84},{-64,-64.84},{-64,-42},{-76,-42}},
            color={127,0,0},
            thickness=0.5));
        connect(respiratoryUnit[i].airways[2], volumeOutflow.q_in) annotation (Line(
            points={{-6,-63.48},{34,-63.48},{34,-42},{64,-42}},
            color={127,0,0},
            thickness=0.5));
        connect(respiratoryUnit[i].blood_out, pulmonaryVeins.q_in[4]) annotation (Line(
            points={{10.32,-86.4},{22,-86.4},{22,-98},{42,-98},{42,-101.025},{41.9,-101.025}},
            color={127,0,0},
            thickness=0.5));
        connect(pressure.y, respiratoryUnit[i].thoraxPressure) annotation (Line(
            points={{-33,-54},{-12.56,-54},{-12.56,-67.68}}, color={0,0,127}));
      end for;
      connect(leftHeartPump.solutionFlow, leftCardiacOutput.y) annotation (
        Line(points = {{49, -150}, {71, -150}}, color = {0, 0, 127}));
      for i in 1:NT loop
        connect(tissueUnit[i].q_in, systemicArteries.q_in[i+1]) annotation (
          Line(points={{13.6611,-195.21},{28,-195.21},{28,-196},{36.1,-196},{
                36.1,-196}},                                                                           color = {127, 0, 0}, thickness = 0.5));
        connect(tissueUnit[i].q_out, systemicVeins.q_in[i+2]) annotation (
          Line(points={{-15.2278,-194.93},{-50.1,-194.93},{-50.1,-196}},        color = {127, 0, 0}, thickness = 0.5));
      end for;
      connect(respiratoryCenter.deadSpaceVentilation, deadSpaceVentilation.solutionFlow)
        annotation (Line(points={{28,2.8},{28,-26},{-4,-26},{-4,-35}},      color={0,
              0,127}));
      connect(respiratoryCenter.ventilation, volumeOutflow.solutionFlow)
        annotation (Line(points={{36.4,2.95},{36.4,-26},{74,-26},{74,-35}},
            color={0,0,127}));

      connect(arterial.pO2, respiratoryCenter.pO2) annotation (Line(points={{58,-165},
              {58,-158},{88,-158},{88,0},{52.9,0},{52.9,5.5}},        color={0,0,127}));
      connect(arterial.pCO2, respiratoryCenter.pCO2) annotation (Line(points={{54,-165},
              {54,-166},{58,-166},{58,-162},{60,-162},{60,-160},{90,-160},{90,17.5},
              {52.6,17.5}},                              color={0,0,127}));
      connect(arterial.b_port, leftHeartPump.q_out) annotation (
        Line(points = {{43.8, -169.8}, {42, -169.8}, {42, -160}}, color = {127, 0, 0}, thickness = 0.5));
      connect(arterial.a_port, systemicArteries.q_in[1]) annotation (
        Line(points={{43.8,-182},{40,-182},{40,-196},{36.1,-196}},              color = {127, 0, 0}, thickness = 0.5));
      connect(systemicVeins.q_in[2], venous.a_port) annotation (
        Line(points={{-50.1,-196},{-50.1,-188},{-55.8,-188},{-55.8,-182}},              color = {127, 0, 0}, thickness = 0.5));
      connect(venous.b_port, rightHeartPump.q_in) annotation (
        Line(points = {{-55.8, -169.8}, {-55.8, -164.9}, {-50, -164.9}, {-50, -160}}, color = {127, 0, 0}, thickness = 0.5));

      connect(systemicArteries.substances, arterial.substances)
        annotation (Line(
          points={{46,-196},{54,-196},{54,-186}},
          color={158,66,200},
          thickness=0.5));
      connect(systemicVeins.substances, venous.substances)
        annotation (Line(
          points={{-60,-196},{-66,-196},{-66,-186}},
          color={158,66,200},
          thickness=0.5));
          //plethy = 120*133.322 * arterial.sO2;//bloodHemodynamics.EithaPressure.pressure * arterial.sO2;
      connect(hemodynamicPressure.y, PPG.u2) annotation (Line(points={{77,-208},{78,
              -208},{78,-204},{76.8,-204},{76.8,-195.6}}, color={0,0,127}));
      connect(PPG.y, plethy) annotation (Line(points={{90.6,-192},{90,-194},{104,-194}},
            color={0,0,127}));
      connect(arterial.sO2, PPG.u1) annotation (Line(points={{62,-165},{62,-162},{74,
              -162},{74,-184},{72,-184},{72,-188.4},{76.8,-188.4}}, color={0,0,127}));
      connect(AirO2.u, AirO2Fraction.y) annotation (Line(points={{-34,20},{-50,20},{
              -50,26},{-53,26}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent={{-180,-220},{100,
                40}}),
            graphics={Rectangle(
              extent={{-86,28},{96,-184}},
              lineColor={28,108,200},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent={{-180,-220},{100,
                40}})),
        experiment(StopTime = 1800, __Dymola_Algorithm = "Dassl"));
    end BloodyMaryPPG2;

    package MeursModel2011
      "models of cardiovascular system used in www.physiome.cz/atlas"
      extends Modelica.Icons.ExamplesPackage;

      package Parts "Utility components used by package KofranekModels2013"
        extends Modelica.Icons.UtilitiesPackage;

        model AtrialElastance
          extends HeartIntervals;
          Physiolibrary.Types.RealIO.HydraulicElastanceOutput Et "elasticity" annotation (
            Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, -20}, {138, 18}})));
          parameter Physiolibrary.Types.HydraulicElastance EMIN "Diastolic elastance";
          parameter Boolean useEs_extInput = false "=true, if external elastance/compliance value is used" annotation (
            Evaluate = true,
            HideResult = true,
            choices(checkBox = true),
            Dialog(group = "Conditional inputs"));
          parameter Physiolibrary.Types.HydraulicElastance EMAX "Maximum systolic elastance" annotation (
            Dialog(enable = not useEs_extInput));
          Physiolibrary.Types.RealIO.HydraulicComplianceInput Es_ext = 1 / es_int if useEs_extInput annotation (
            Placement(transformation(extent = {{60, 60}, {100, 100}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, origin = {-80, 80})));
        protected
          Physiolibrary.Types.HydraulicElastance es_int;
        equation
          if not useEs_extInput then
            es_int = EMAX;
          end if;
          Et = smooth(1, if time - T0 < Tas then EMIN + (es_int - EMIN) * sin(Modelica.Constants.pi * (time - T0) / Tas) else EMIN);
          annotation (
            Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-100, 82}, {100, -100}}, pattern = LinePattern.None, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, lineColor = {0, 0, 255}), Text(extent = {{-98, 82}, {98, 24}}, lineColor = {0, 0, 255}, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "Atrial elastance"), Line(points = {{-78, -34}, {-76, -26}, {-70, -14}, {-58, 6}, {-36, 36}, {-14, 14}, {-6, -10}, {0, -32}, {6, -34}, {88, -34}, {94, -34}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Text(extent = {{-220, -102}, {200, -120}}, lineColor = {0, 0, 255}, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "%name"), Text(extent = {{72, 4}, {102, -8}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "Ct")}));
        end AtrialElastance;

        model VentricularElastance
          extends HeartIntervals;
          Physiolibrary.Types.RealIO.HydraulicElastanceOutput Et "ventricular elasticity" annotation (
            Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, 4}, {138, 42}})));
          Modelica.Blocks.Interfaces.RealOutput Et0 "normalized ventricular elasticity (0..1)" annotation (
            Placement(transformation(extent = {{100, -24}, {120, -4}}), iconTransformation(extent = {{100, -40}, {138, -2}})));
          Physiolibrary.Types.RealIO.TimeOutput HeartInterval "eapsed time" annotation (
            Placement(transformation(extent = {{102, -42}, {122, -22}}), iconTransformation(extent = {{100, -98}, {138, -60}})));
          parameter Physiolibrary.Types.HydraulicElastance EMIN "Diastolic elastance ";
          constant Real Kn = 0.57923032735652;
          //Kn is always = 0.5792303273565197
          //... the t * sin(pi*t) has its maximum at t = 0.645773676543406 and = 0.5792303273565197
          //Equation to calculate normalized elastance ET0 was:
          //Et0=EMIN+(EMAX-EMIN)*((time-T0)-(Tas+Tav))/Tvs)*sin(Modelica.Constants.pi*(((time-T0)-(Tas+Tav))/Tvs));
          parameter Boolean useEs_extInput = false "=true, if external elastance/compliance value is used" annotation (
            Evaluate = true,
            HideResult = true,
            choices(checkBox = true),
            Dialog(group = "Conditional inputs"));
          parameter Physiolibrary.Types.HydraulicElastance EMAX "Maximum systolic elastance" annotation (
            Dialog(enable = not useEs_extInput));
          Physiolibrary.Types.RealIO.HydraulicComplianceInput Es_ext(start = 1 / EMAX) = 1 / es_int if useEs_extInput annotation (
            Placement(transformation(extent = {{60, 60}, {100, 100}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, origin = {-80, 80})));
        protected
          Physiolibrary.Types.HydraulicElastance es_int;
        equation
          if not useEs_extInput then
            es_int = EMAX;
          end if;
          HeartInterval = time - T0;
          Et = EMIN + (es_int - EMIN) * Et0;
          Et0 = smooth(1, if HeartInterval >= Tas + Tav and HeartInterval < Tas + Tav + Tvs then (HeartInterval - (Tas + Tav)) / Tvs * sin(Modelica.Constants.pi * (HeartInterval - (Tas + Tav)) / Tvs) / Kn else 0);
          annotation (
            Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-98, 82}, {100, -100}}, pattern = LinePattern.None, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, lineColor = {0, 0, 255}), Text(extent = {{-82, 82}, {80, 24}}, lineColor = {0, 0, 255}, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "Ventricular elastance"), Line(points = {{-72, -34}, {-62, -34}, {-52, -34}, {-44, 8}, {-18, 38}, {-12, 14}, {-6, -10}, {0, -32}, {6, -34}, {88, -34}, {94, -34}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Text(extent = {{-220, -102}, {200, -120}}, lineColor = {0, 0, 255}, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "%name"), Text(extent = {{96, -32}, {68, -8}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "Et0"), Text(extent = {{42, -72}, {88, -84}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "Heart interval"), Text(extent = {{62, 30}, {96, 8}}, lineColor = {0, 0, 255}, lineThickness = 1, fillColor = {255, 255, 170}, fillPattern = FillPattern.Solid, textString = "Ct")}));
        end VentricularElastance;

        model HeartIntervals
          discrete Physiolibrary.Types.Time Tas;
          discrete Physiolibrary.Types.Time T0;
          discrete Physiolibrary.Types.Time Tvs;
          parameter Physiolibrary.Types.Time Tav(displayUnit = "s") = 0.01 "atrioventricular delay";
          parameter Physiolibrary.Types.Frequency HR_start = 1.2 "initial heart rate";
          discrete Modelica.Units.SI.Time HP(start = 0) "heart period";
          Boolean b(start = false);
          Physiolibrary.Types.RealIO.FrequencyInput HR(start = 1.2) "heart rate" annotation (
            Placement(transformation(extent = {{-12, 68}, {28, 108}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 80})));
        initial equation
          T0 = 0 "start time of cardiac cycle";
          HP = 1 / HR_start "update heart period per heart rate";
          Tas = 0.03 + 0.09 / HR_start "duration of atrial systole";
          Tvs = 0.16 + 0.2 / HR_start "duration of ventricular systole";
        equation
          b = time - pre(T0) >= pre(HP) "true if new pulse occurs";
          when {b} then
            T0 = time "start time of cardiac cycle";
            HP = 1 / HR "update heart period per heart rate";
            Tas = 0.03 + 0.09 * HP "duration of atrial systole";
            Tvs = 0.16 + 0.2 * HP "duration of ventricular systole";
          end when;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-64, 102}, {-6, 78}}, lineColor = {0, 0, 255}, textString = "HR")}));
        end HeartIntervals;

        model ECGGenerator
          Real ecgwave[5]; //ecg curve parts separate P,Q,R,S,T waves
          Real dt[5];
          Real ecg; //ecg curve, sum of PQRST waves
          Real relt; //relative time
          discrete Physiolibrary.Types.Time T0; //begining time of next heartinterval
          parameter Physiolibrary.Types.Frequency HR_start = 1.2 "initial heart rate";
          discrete Modelica.Units.SI.Time HP(start = 0) "heart period";
          Boolean b(start = false);
          Physiolibrary.Types.RealIO.FrequencyInput HR(start = 1.2) "heart rate" annotation (
            Placement(transformation(extent = {{-12, 68}, {28, 108}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 80})));
          parameter Real A[5] = {3.6,-69.5,142,-105,4};
          parameter Real B[5] = {0.24,0.1,0.11,0.11,0.42};
          parameter Real d[5] = {-69.5+360,-12.5,0,28.5,146}; //P wave +360 to match hemodynamics heartintervals
          parameter Real drad[5] =  { (d[i]+90)/180*Modelica.Constants.pi for i in 1:5};
        initial equation
          T0 = 0 "start time of cardiac cycle";
          HP = 1 / HR_start "update heart period per heart rate";
        equation
          b = time - pre(T0) >= pre(HP) "true if new pulse occurs";
          when {b} then
            T0 = time "start time of cardiac cycle";
            HP = 1 / HR "update heart period per heart rate";

          end when;

          relt = (time-T0)/HP;  //relative time in period

          for i in 1:5 loop
            dt[i] =  relt*2*Modelica.Constants.pi - drad[i] + Modelica.Constants.pi/4; //+pi/4 to match hemodynamics QRS just before ventricular systole
            ecgwave[i] = A[i]*dt[i] *
                         exp( - dt[i]^2 / (2*B[i]^2));
          end for;
          der(ecg) = if (HP>20) then 0
          else - sum(ecgwave) - ecg;
          annotation (Icon(graphics={Rectangle(
                  extent={{-82,62},{80,-82}},
                  lineColor={28,108,200},
                  fillColor={255,255,170},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-110,-80},{138,-106}},
                  textColor={28,108,200},
                  textString="%name"),
                Polygon(
                  points={{-34,-6},{-24,-6},{-22,-2},{-20,-6},{-16,-10},{-12,26},{-10,-36},
                      {-8,0},{-6,-10},{-2,-6},{6,-6},{16,2},{22,2},{32,-6},{-34,-6}},
                  lineColor={28,108,200},
                  fillColor={255,255,170},
                  fillPattern=FillPattern.Solid)}));
        end ECGGenerator;

        model TestEcg
          ECGGenerator eCGGenerator
            annotation (Placement(transformation(extent={{8,-30},{28,-10}})));
          replaceable Physiolibrary.Types.Constants.FrequencyConst HeartRate(k(
                displayUnit="1/min") = 1.2)                                                                  annotation (
            Placement(transformation(origin={-49,4.5},       extent = {{-11, -6.5}, {11, 6.5}})));
        equation
          connect(HeartRate.y, eCGGenerator.HR) annotation (Line(points={{
                  -35.25,4.5},{18,4.5},{18,-12}}, color={0,0,127}));
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=false)),
            Diagram(coordinateSystem(preserveAspectRatio=false)),
            experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=5000,
              Tolerance=1e-05,
              __Dymola_Algorithm="Dassl"));
        end TestEcg;
      end Parts;

      model Hemodynamics
        extends Physiolibrary.Icons.CardioVascular;
        Physiolibrary.Fluid.Components.ElasticVesselElastance Epa(
        volume_start=0.000106,
        ZeroPressureVolume=5e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,
        Elastance=31064116.267695,                                                                                                                                                            nPorts = 2) annotation (
          Placement(transformation(extent = {{-94, 84}, {-66, 112}})));
        Physiolibrary.Fluid.Components.Resistor Rpp(Resistance(displayUnit = "(mmHg.s)/ml") = 14665462.61565) annotation (
          Placement(transformation(extent = {{-56, 85}, {-22, 111}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Epv(
        volume_start=0.000518,
        ZeroPressureVolume=0.00035,
        ExternalPressure=-533.28954966,
        Elastance=6066168.6273825,                                                                                                                                                              nPorts = 2) annotation (
          Placement(transformation(extent = {{-14, 84}, {20, 112}})));
        Physiolibrary.Fluid.Components.Resistor Rlain(Resistance(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{26, 86}, {56, 110}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance LeftAtrium(
        isExternalPressureAbsolute=false,                                useElastanceInput = true,
        volume_start=9.31e-05,
        ZeroPressureVolume=3e-05,
        ExternalPressure=-533.28954966,                                                                                                                                                           nPorts = 2) annotation (
          Placement(transformation(extent = {{74, 50}, {102, 78}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance LeftVentricle(
        isExternalPressureAbsolute=false,                                   useElastanceInput = true,
        volume_start=0.000144,
        ZeroPressureVolume=6e-05,
        ExternalPressure=-533.28954966,                                                                                                                                                              nPorts = 2) annotation (
          Placement(transformation(extent = {{150, 50}, {178, 78}})));
        Physiolibrary.Fluid.Components.IdealValveResistance AorticValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 1066579.09932) annotation (
          Placement(transformation(extent = {{184, 76}, {208, 52}})));
        Parts.AtrialElastance LAtrialElastance(Tav(displayUnit = "s"), EMIN = 15998686.4898, EMAX = 37330268.4762) annotation (
          Placement(transformation(extent = {{80, 92}, {118, 124}})));
        Parts.VentricularElastance LVentricularElastance(EMIN = 11999014.86735, EMAX = 533289549.66) annotation (
          Placement(transformation(extent = {{164, 88}, {200, 120}})));
        Physiolibrary.Fluid.Components.IdealValveResistance MitralValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(origin = {127, 64}, extent = {{-13, 12}, {13, -12}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eitha(
        ExternalPressure=-533.28954966,                                                                       nPorts = 4,
        volume_start=0.0002,
        ZeroPressureVolume=0.00014,
        Elastance=190651014.00345)                                                                                                                                                                          annotation (
          Placement(transformation(extent = {{168, 6}, {190, 28}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eetha(
        volume_start=0.000526,
        ZeroPressureVolume=0.00037,
        Elastance=74127247.40274,                                                                                                                      nPorts = 3) annotation (
          Placement(transformation(extent = {{56, 4}, {82, 30}})));
        Physiolibrary.Fluid.Components.Inertia inertia(I(displayUnit=
              "mmHg.s2/g") = 226.6480586055, massFlow_start(displayUnit=
              "g/min") = 0.021666666666667)                                                                                                                       annotation (
          Placement(transformation(extent = {{-11, -11}, {11, 11}}, rotation = 180, origin = {141, 17})));
        Physiolibrary.Fluid.Components.Resistor Retha(Resistance(displayUnit = "(mmHg.s)/ml") = 7999343.2449) annotation (
          Placement(transformation(extent = {{90, 6}, {112, 28}})));
        Physiolibrary.Fluid.Components.Resistor Rsart(Resistance(displayUnit = "(mmHg.s)/ml") = 106657909.932) annotation (
          Placement(transformation(extent = {{14, -13}, {-14, 13}}, origin = {24, 19})));
        Physiolibrary.Fluid.Components.Resistor Rsven(Resistance(displayUnit = "(mmHg.s)/ml") = 26664477.483) annotation (
          Placement(transformation(extent = {{14, -13}, {-14, 13}}, origin = {-60, 17})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Est(volume_start = 0.283e-3, ZeroPressureVolume = 0.185e-3, Elastance = 34930465.50273, nPorts = 3) annotation (
          Placement(transformation(extent = {{-28, 6}, {-4, 28}})));
        Physiolibrary.Fluid.Components.Resistor Rethv(Resistance(displayUnit = "(mmHg.s)/ml") = 11999014.86735) annotation (
          Placement(transformation(extent = {{-120, 4}, {-146, 30}})));
        Physiolibrary.Fluid.Components.Resistor Rrain(Resistance(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{-208, 4}, {-236, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eithv(
        volume_start=0.00148,
        ZeroPressureVolume=0.00119,
        ExternalPressure=-533.28954966,
        Elastance=2426467.450953,                                                                                                                                                               nPorts = 3) annotation (
          Placement(transformation(extent = {{-194, 4}, {-166, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eethv(
        volume_start=0.00153,
        ZeroPressureVolume=0.001,
        Elastance=2253148.3473135,                                                                                                                  nPorts = 3) annotation (
          Placement(transformation(extent = {{-108, 4}, {-82, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance RightAtrium(
        isExternalPressureAbsolute=false,                                 useElastanceInput = true,
        volume_start=0.000135,
        ZeroPressureVolume=3e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,                                                                                                                                      nPorts = 2) annotation (
          Placement(transformation(extent = {{-242, 44}, {-214, 72}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance RightVentricle(
        isExternalPressureAbsolute=false,                                    useElastanceInput = true,
        volume_start=0.000131,
        ZeroPressureVolume=4e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,                                                                                                                                         nPorts = 2) annotation (
          Placement(transformation(extent = {{-170, 42}, {-140, 72}})));
        Physiolibrary.Fluid.Components.IdealValveResistance PulmonaryValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{-132, 70}, {-106, 44}})));
        Parts.AtrialElastance RAtrialElastance(EMIN = 6666119.37075, EMAX = 19998358.11225) annotation (
          Placement(transformation(extent = {{-244, 86}, {-206, 118}})));
        Parts.VentricularElastance RVentricularElastance(EMIN = 7599376.082655, EMAX = 65327969.83335) annotation (
          Placement(transformation(extent = {{-180, 88}, {-150, 122}})));
        Physiolibrary.Fluid.Components.IdealValveResistance TricuspidValve(_Roff = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(origin = {-189, 58}, extent = {{-13, 12}, {13, -12}})));
        inner Modelica.Fluid.System system(p_ambient(displayUnit = "mmHg") = 101325.0144354, T_ambient = 310.15) annotation (
          Placement(transformation(extent = {{-48, 144}, {-28, 164}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EithaPressure
          annotation (Placement(transformation(extent={{190,-36},{210,-16}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EethaPressure
          annotation (Placement(transformation(extent={{78,-56},{98,-36}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EstPressure
          annotation (Placement(transformation(extent={{-12,-52},{8,-32}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EethvPressure
          annotation (Placement(transformation(extent={{-92,-52},{-72,-32}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EithvPressure annotation (
            Placement(transformation(extent={{-166,-54},{-146,-34}})));
        Physiolibrary.Fluid.Sensors.Power power annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-96,72})));
        Physiolibrary.Fluid.Sensors.Power power1 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-254,36})));
        Modelica.Blocks.Math.Mean rightHeartPower(f(displayUnit = "1/min") = 1.2) annotation (
          Placement(transformation(extent = {{-214, -66}, {-194, -46}})));
        Modelica.Blocks.Math.Feedback feedback annotation (
          Placement(transformation(extent = {{-264, -46}, {-244, -66}})));
        Physiolibrary.Fluid.Sensors.Power power2 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,82})));
        Physiolibrary.Fluid.Sensors.Power power3 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={224,38})));
        Modelica.Blocks.Math.Mean leftHeartPower(f(displayUnit = "1/min") = 1.2) annotation (
          Placement(transformation(extent = {{186, -110}, {206, -90}})));
        Modelica.Blocks.Math.Feedback feedback1 annotation (
          Placement(transformation(extent = {{136, -110}, {156, -90}})));
        Physiolibrary.Fluid.Sensors.Sphygmomanometer arterialPressure(
            MeasurementTime(displayUnit="s") = 60/72)
          annotation (Placement(transformation(extent={{198,-66},{218,-46}})));
        Modelica.Blocks.Noise.UniformNoise uniformNoise(
          samplePeriod=6,
          y_min=-0.06,
          y_max=0.06) annotation (Placement(transformation(extent={{-274,164},{
                  -254,184}})));
        Modelica.Blocks.Math.Add randomizedHeartRate annotation (Placement(
              transformation(extent={{-236,158},{-216,178}})));
        Modelica.Blocks.Logical.Switch currentHeartReat annotation (Placement(
              transformation(extent={{-182,138},{-162,158}})));
        Modelica.Blocks.Sources.BooleanConstant randomizeHR annotation (
            Placement(transformation(extent={{-214,142},{-194,162}})));
        Parts.ECGGenerator Ecg
          annotation (Placement(transformation(extent={{-100,146},{-80,166}})));
        Physiolibrary.Types.Constants.FrequencyConst heartRate(k=1.2)
          annotation (Placement(transformation(extent={{-274,132},{-266,140}})));
      equation
        connect(Retha.q_out, inertia.q_out) annotation (
          Line(points = {{112, 17}, {130, 17}}, thickness = 1));
        connect(RightAtrium.q_in[1], TricuspidValve.q_in) annotation (
          Line(points={{-228.14,57.09},{-214,57.09},{-214,58},{-202,58}},          color = {127, 0, 0}, thickness = 0.5));
        connect(TricuspidValve.q_out, RightVentricle.q_in[1]) annotation (
          Line(points={{-176,58},{-166,58},{-166,56.025},{-155.15,56.025}},        color = {127, 0, 0}, thickness = 0.5));
        connect(RightVentricle.q_in[2], PulmonaryValve.q_in) annotation (
          Line(points={{-155.15,57.975},{-143.725,57.975},{-143.725,57},{-132,
              57}},                                                                        color = {127, 0, 0}, thickness = 0.5));
        connect(Epa.q_in[1], Rpp.q_in) annotation (
          Line(points={{-80.14,97.09},{-68,97.09},{-68,98},{-56,98}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rpp.q_out, Epv.q_in[1]) annotation (
          Line(points={{-22,98},{-8,98},{-8,97.09},{2.83,97.09}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Epv.q_in[2], Rlain.q_in) annotation (
          Line(points={{2.83,98.91},{16,98.91},{16,98},{26,98}},          color = {127, 0, 0}, thickness = 0.5));
        connect(LeftAtrium.q_in[1], MitralValve.q_in) annotation (
          Line(points={{87.86,63.09},{102,63.09},{102,64},{114,64}},          color = {127, 0, 0}, thickness = 0.5));
        connect(MitralValve.q_out, LeftVentricle.q_in[1]) annotation (
          Line(points={{140,64},{154,64},{154,66},{163.86,66},{163.86,63.09}},            color = {127, 0, 0}, thickness = 0.5));
        connect(LeftVentricle.q_in[2], AorticValve.q_in) annotation (
          Line(points={{163.86,64.91},{172,64.91},{172,64},{184,64}},          color = {127, 0, 0}, thickness = 0.5));
        connect(inertia.q_in, Eitha.q_in[1]) annotation (
          Line(points={{152,17},{164,17},{164,15.9275},{178.89,15.9275}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Retha.q_in, Eetha.q_in[1]) annotation (
          Line(points={{90,17},{80,17},{80,15.8733},{68.87,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rsart.q_in, Eetha.q_in[2]) annotation (
          Line(points = {{38, 19}, {52, 19}, {52, 17}, {68.87, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Est.q_in[1], Rsart.q_out) annotation (
          Line(points={{-16.12,16.0467},{-3.18,16.0467},{-3.18,19},{10,19}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Est.q_in[2], Rsven.q_in) annotation (
          Line(points = {{-16.12, 17}, {-31.18, 17}, {-31.18, 17}, {-46, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Rsven.q_out, Eethv.q_in[1]) annotation (
          Line(points={{-74,17},{-86,17},{-86,15.8733},{-95.13,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rethv.q_in, Eethv.q_in[2]) annotation (
          Line(points = {{-120, 17}, {-107, 17}, {-107, 17}, {-95.13, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Rethv.q_out, Eithv.q_in[1]) annotation (
          Line(points={{-146,17},{-164,17},{-164,15.8733},{-180.14,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rrain.q_in, Eithv.q_in[2]) annotation (
          Line(points = {{-208, 17}, {-194, 17}, {-194, 17}, {-180.14, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(EithaPressure.port, Eitha.q_in[2]) annotation (
          Line(points={{200,-36},{178.89,-36},{178.89,16.6425}},       color = {127, 0, 0}, thickness = 0.5));
        connect(EethaPressure.port, Eetha.q_in[3]) annotation (
          Line(points={{88,-56},{88,-54},{68.87,-54},{68.87,18.1267}},          color = {127, 0, 0}, thickness = 0.5));
        connect(EstPressure.port, Est.q_in[3]) annotation (
          Line(points={{-2,-52},{-16.12,-52},{-16.12,17.9533}},        color = {127, 0, 0}, thickness = 0.5));
        connect(EethvPressure.port, Eethv.q_in[3]) annotation (
          Line(points={{-82,-52},{-95.13,-52},{-95.13,18.1267}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Eithv.q_in[3],EithvPressure.port)  annotation (
          Line(points={{-180.14,18.1267},{-180.14,-54},{-156,-54}},        color = {127, 0, 0}, thickness = 0.5));
        connect(RAtrialElastance.Et, RightAtrium.elastance) annotation (
          Line(points = {{-202.39, 101.84}, {-202.39, 85.92}, {-225.2, 85.92}, {-225.2, 70.6}}, color = {0, 0, 127}));
        connect(RVentricularElastance.Et, RightVentricle.elastance) annotation (
          Line(points = {{-147.15, 108.91}, {-134, 108.91}, {-134, 76}, {-152, 76}, {-152, 70.5}}, color = {0, 0, 127}));
        connect(LAtrialElastance.Et, LeftAtrium.elastance) annotation (
          Line(points = {{121.61, 107.84}, {121.61, 91.92}, {90.8, 91.92}, {90.8, 76.6}}, color = {0, 0, 127}));
        connect(LVentricularElastance.Et, LeftVentricle.elastance) annotation (
          Line(points = {{203.42, 107.68}, {222, 107.68}, {222, 76.6}, {166.8, 76.6}}, color = {0, 0, 127}));
        connect(PulmonaryValve.q_out, power.q_in) annotation (
          Line(points = {{-106, 57}, {-96, 57}, {-96, 62}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power.q_out, Epa.q_in[2]) annotation (
          Line(points={{-96,82},{-96,98.91},{-80.14,98.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Rrain.q_out, power1.q_in) annotation (
          Line(points = {{-236, 17}, {-254, 17}, {-254, 26}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power1.q_out, RightAtrium.q_in[2]) annotation (
          Line(points={{-254,46},{-254,58.91},{-228.14,58.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(power.power, feedback.u1) annotation (
          Line(points = {{-108, 72}, {-276, 72}, {-276, -56}, {-262, -56}}, color = {0, 0, 127}));
        connect(power1.power, feedback.u2) annotation (
          Line(points = {{-266, 36}, {-268, 36}, {-268, -30}, {-254, -30}, {-254, -48}}, color = {0, 0, 127}));
        connect(feedback.y, rightHeartPower.u) annotation (
          Line(points = {{-245, -56}, {-216, -56}}, color = {0, 0, 127}));
        connect(Rlain.q_out, power2.q_in) annotation (
          Line(points = {{56, 98}, {70, 98}, {70, 92}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power2.q_out, LeftAtrium.q_in[2]) annotation (
          Line(points={{70,72},{70,64.91},{87.86,64.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(AorticValve.q_out, power3.q_in) annotation (
          Line(points = {{208, 64}, {224, 64}, {224, 48}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power3.q_out, Eitha.q_in[3]) annotation (
          Line(points={{224,28},{224,17.3575},{178.89,17.3575}},      color = {127, 0, 0}, thickness = 0.5));
        connect(feedback1.y, leftHeartPower.u) annotation (
          Line(points = {{155, -100}, {184, -100}}, color = {0, 0, 127}));
        connect(power3.power, feedback1.u1) annotation (
          Line(points = {{236, 38}, {244, 38}, {244, -84}, {130, -84}, {130, -100}, {138, -100}}, color = {0, 0, 127}));
        connect(power2.power, feedback1.u2) annotation (
          Line(points = {{82, 82}, {106, 82}, {106, -8}, {118, -8}, {118, -124}, {146, -124}, {146, -108}}, color = {0, 0, 127}));
        connect(Eitha.q_in[4],arterialPressure.port)  annotation (
          Line(points={{178.89,18.0725},{178.89,-66},{208,-66}},       color = {127, 0, 0}, thickness = 0.5));
        connect(uniformNoise.y, randomizedHeartRate.u1)
          annotation (Line(points={{-253,174},{-238,174}}, color={0,0,127}));
        connect(randomizedHeartRate.y, currentHeartReat.u1) annotation (Line(
              points={{-215,168},{-184,168},{-184,156}}, color={0,0,127}));
        connect(currentHeartReat.y, RVentricularElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,134},{-165,134},{-165,118.6}},
              color={0,0,127}));
        connect(currentHeartReat.y, LAtrialElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,136},{99,136},{99,120.8}},
              color={0,0,127}));
        connect(currentHeartReat.y, RAtrialElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,130},{-225,130},{-225,114.8}},
              color={0,0,127}));
        connect(currentHeartReat.y, LVentricularElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,136},{182,136},{182,116.8}},
              color={0,0,127}));
        connect(randomizeHR.y, currentHeartReat.u2) annotation (Line(points={{
                -193,152},{-190,152},{-190,148},{-184,148}}, color={255,0,255}));
        connect(Ecg.HR, currentHeartReat.y) annotation (Line(points={{-90,164},
                {-90,172},{-152,172},{-152,148},{-161,148}}, color={0,0,127}));
        connect(heartRate.y, currentHeartReat.u3) annotation (Line(points={{
                -265,136},{-192,136},{-192,132},{-184,132},{-184,140}}, color={
                0,0,127}));
        connect(heartRate.y, randomizedHeartRate.u2) annotation (Line(points={{
                -265,136},{-248,136},{-248,162},{-238,162}}, color={0,0,127}));
        annotation (
          Diagram(coordinateSystem(extent = {{-280, -140}, {280, 180}}, preserveAspectRatio = false)),
          Icon(coordinateSystem(extent = {{-280, -140}, {280, 180}}, preserveAspectRatio = false), graphics),
          Documentation(info = "<html>
        <p>Model of cardiovascular system using to demonstrate elastic and resistance features of veins and arteries in pulmonary and systemic circulation and influence of cardiac output on it.</p>
        <ul>
        <li>J. A. Goodwin, W. L. van Meurs, C. D. Sa Couto, J. E. W.Beneken, S. A. Graves, A model for educational simulation of infant cardiovascular physiology., Anesthesia and analgesia 99 (6)(2004) 1655&ndash;1664. doi:10.1213/01.ANE.0000134797.52793.AF.</li>
        <li>C. D. Sa Couto, W. L. van Meurs, J. A. Goodwin, P. Andriessen,A Model for Educational Simulation of Neonatal Cardiovascular Pathophysiology, Simulation in Healthcare 1 (Inaugural) (2006) 4&ndash;12.</li>
        <li>W. van Meurs, Modeling and Simulation in Biomedical Engineering: Applications in Cardiorespiratory Physiology, McGraw-Hill Professional, 2011.</li>
        </ul>
        </html>",
                revisions = "<html>
        <ul>
        <li><i>Jul 2015 </i>by Tomas Kulhanek: Created. </li>
        </ul>
        </html>"),
          experiment(StopTime = 10, __Dymola_Algorithm = "Dassl"));
      end Hemodynamics;

      model testHemodynamics
        extends Physiolibrary.Icons.CardioVascular;
        Physiolibrary.Fluid.Components.ElasticVesselElastance Epa(
        volume_start=0.000106,
        ZeroPressureVolume=5e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,
        Elastance=31064116.267695,                                                                                                                                                            nPorts = 2) annotation (
          Placement(transformation(extent = {{-94, 84}, {-66, 112}})));
        Physiolibrary.Fluid.Components.Resistor Rpp(Resistance(displayUnit = "(mmHg.s)/ml") = 14665462.61565) annotation (
          Placement(transformation(extent = {{-56, 85}, {-22, 111}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Epv(
        volume_start=0.000518,
        ZeroPressureVolume=0.00035,
        ExternalPressure=-533.28954966,
        Elastance=6066168.6273825,                                                                                                                                                              nPorts = 2) annotation (
          Placement(transformation(extent = {{-14, 84}, {20, 112}})));
        Physiolibrary.Fluid.Components.Resistor Rlain(Resistance(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{26, 86}, {56, 110}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance LeftAtrium(
        isExternalPressureAbsolute=false,                                useElastanceInput = true,
        volume_start=9.31e-05,
        ZeroPressureVolume=3e-05,
        ExternalPressure=-533.28954966,                                                                                                                                                           nPorts = 2) annotation (
          Placement(transformation(extent = {{74, 50}, {102, 78}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance LeftVentricle(
        isExternalPressureAbsolute=false,                                   useElastanceInput = true,
        volume_start=0.000144,
        ZeroPressureVolume=6e-05,
        ExternalPressure=-533.28954966,                                                                                                                                                              nPorts = 2) annotation (
          Placement(transformation(extent = {{150, 50}, {178, 78}})));
        Physiolibrary.Fluid.Components.IdealValveResistance AorticValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 1066579.09932) annotation (
          Placement(transformation(extent = {{184, 76}, {208, 52}})));
        Parts.AtrialElastance LAtrialElastance(Tav(displayUnit = "s"), EMIN = 15998686.4898, EMAX = 37330268.4762) annotation (
          Placement(transformation(extent = {{80, 92}, {118, 124}})));
        Parts.VentricularElastance LVentricularElastance(EMIN = 11999014.86735, EMAX = 533289549.66) annotation (
          Placement(transformation(extent = {{164, 88}, {200, 120}})));
        Physiolibrary.Fluid.Components.IdealValveResistance MitralValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(origin = {127, 64}, extent = {{-13, 12}, {13, -12}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eitha(
        ExternalPressure=-533.28954966,                                                                       nPorts = 4,
        volume_start=0.0002,
        ZeroPressureVolume=0.00014,
        Elastance=190651014.00345)                                                                                                                                                                          annotation (
          Placement(transformation(extent = {{168, 6}, {190, 28}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eetha(
        volume_start=0.000526,
        ZeroPressureVolume=0.00037,
        Elastance=74127247.40274,                                                                                                                      nPorts = 3) annotation (
          Placement(transformation(extent = {{56, 4}, {82, 30}})));
        Physiolibrary.Fluid.Components.Inertia inertia(I(displayUnit=
              "mmHg.s2/g") = 226.6480586055, massFlow_start(displayUnit=
              "g/min") = 0.021666666666667)                                                                                                                       annotation (
          Placement(transformation(extent = {{-11, -11}, {11, 11}}, rotation = 180, origin = {141, 17})));
        Physiolibrary.Fluid.Components.Resistor Retha(Resistance(displayUnit = "(mmHg.s)/ml") = 7999343.2449) annotation (
          Placement(transformation(extent = {{90, 6}, {112, 28}})));
        Physiolibrary.Fluid.Components.Resistor Rsart(Resistance(displayUnit = "(mmHg.s)/ml") = 106657909.932) annotation (
          Placement(transformation(extent = {{14, -13}, {-14, 13}}, origin = {24, 19})));
        Physiolibrary.Fluid.Components.Resistor Rsven(Resistance(displayUnit = "(mmHg.s)/ml") = 26664477.483) annotation (
          Placement(transformation(extent = {{14, -13}, {-14, 13}}, origin = {-60, 17})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Est(volume_start = 0.283e-3, ZeroPressureVolume = 0.185e-3, Elastance = 34930465.50273, nPorts = 3) annotation (
          Placement(transformation(extent = {{-28, 6}, {-4, 28}})));
        Physiolibrary.Fluid.Components.Resistor Rethv(Resistance(displayUnit = "(mmHg.s)/ml") = 11999014.86735) annotation (
          Placement(transformation(extent = {{-120, 4}, {-146, 30}})));
        Physiolibrary.Fluid.Components.Resistor Rrain(Resistance(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{-208, 4}, {-236, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eithv(
        volume_start=0.00148,
        ZeroPressureVolume=0.00119,
        ExternalPressure=-533.28954966,
        Elastance=2426467.450953,                                                                                                                                                               nPorts = 3) annotation (
          Placement(transformation(extent = {{-194, 4}, {-166, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eethv(
        volume_start=0.00153,
        ZeroPressureVolume=0.001,
        Elastance=2253148.3473135,                                                                                                                  nPorts = 3) annotation (
          Placement(transformation(extent = {{-108, 4}, {-82, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance RightAtrium(
        isExternalPressureAbsolute=false,                                 useElastanceInput = true,
        volume_start=0.000135,
        ZeroPressureVolume=3e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,                                                                                                                                      nPorts = 2) annotation (
          Placement(transformation(extent = {{-242, 44}, {-214, 72}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance RightVentricle(
        isExternalPressureAbsolute=false,                                    useElastanceInput = true,
        volume_start=0.000131,
        ZeroPressureVolume=4e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,                                                                                                                                         nPorts = 2) annotation (
          Placement(transformation(extent = {{-170, 42}, {-140, 72}})));
        Physiolibrary.Fluid.Components.IdealValveResistance PulmonaryValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{-132, 70}, {-106, 44}})));
        Parts.AtrialElastance RAtrialElastance(EMIN = 6666119.37075, EMAX = 19998358.11225) annotation (
          Placement(transformation(extent = {{-244, 86}, {-206, 118}})));
        Parts.VentricularElastance RVentricularElastance(EMIN = 7599376.082655, EMAX = 65327969.83335) annotation (
          Placement(transformation(extent = {{-180, 88}, {-150, 122}})));
        Physiolibrary.Fluid.Components.IdealValveResistance TricuspidValve(_Roff = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(origin = {-189, 58}, extent = {{-13, 12}, {13, -12}})));
        inner Modelica.Fluid.System system(p_ambient(displayUnit = "mmHg") = 101325.0144354, T_ambient = 310.15) annotation (
          Placement(transformation(extent = {{-48, 144}, {-28, 164}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EithaPressure
          annotation (Placement(transformation(extent={{190,-36},{210,-16}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EethaPressure
          annotation (Placement(transformation(extent={{78,-56},{98,-36}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EstPressure
          annotation (Placement(transformation(extent={{-12,-52},{8,-32}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EethvPressure
          annotation (Placement(transformation(extent={{-92,-52},{-72,-32}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EithvPressure annotation (
            Placement(transformation(extent={{-166,-54},{-146,-34}})));
        Physiolibrary.Fluid.Sensors.Power power annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-96,72})));
        Physiolibrary.Fluid.Sensors.Power power1 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-254,36})));
        Modelica.Blocks.Math.Mean rightHeartPower(f(displayUnit = "1/min") = 1.2) annotation (
          Placement(transformation(extent = {{-214, -66}, {-194, -46}})));
        Modelica.Blocks.Math.Feedback feedback annotation (
          Placement(transformation(extent = {{-264, -46}, {-244, -66}})));
        Physiolibrary.Fluid.Sensors.Power power2 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,82})));
        Physiolibrary.Fluid.Sensors.Power power3 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={224,38})));
        Modelica.Blocks.Math.Mean leftHeartPower(f(displayUnit = "1/min") = 1.2) annotation (
          Placement(transformation(extent = {{186, -110}, {206, -90}})));
        Modelica.Blocks.Math.Feedback feedback1 annotation (
          Placement(transformation(extent = {{136, -110}, {156, -90}})));
        Physiolibrary.Fluid.Sensors.Sphygmomanometer arterialPressure(
            MeasurementTime(displayUnit="s") = 60/72)
          annotation (Placement(transformation(extent={{198,-66},{218,-46}})));
        Modelica.Blocks.Noise.UniformNoise uniformNoise(
          samplePeriod=6,
          y_min=-0.03,
          y_max=0.03) annotation (Placement(transformation(extent={{-274,164},{
                  -254,184}})));
        Modelica.Blocks.Math.Add randomizedHeartRate annotation (Placement(
              transformation(extent={{-236,158},{-216,178}})));
        Modelica.Blocks.Logical.Switch currentHeartReat annotation (Placement(
              transformation(extent={{-182,138},{-162,158}})));
        Modelica.Blocks.Sources.BooleanConstant randomizeHR annotation (
            Placement(transformation(extent={{-214,142},{-194,162}})));
        Parts.ECGGenerator Ecg
          annotation (Placement(transformation(extent={{-100,146},{-80,166}})));
        Modelica.Blocks.Tables.CombiTable1Ds convertSaturationToHeartRate(table
            =[0,0.03; 0.55,0.03; 0.6,180/60; 0.7,120/60; 0.8,95/60; 0.9,72/60],
            extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
          annotation (Placement(transformation(extent={{-276,120},{-256,140}})));
        Modelica.Blocks.Sources.TimeTable saturationInTime(table=[0,0.9; 10,
              0.90; 20,0.82; 30,0.75; 40,0.63; 50,0.6; 60,0.55; 120,0.55; 130,
              0.65; 140,0.75; 150,0.75]) annotation (Placement(transformation(
                extent={{-312,118},{-292,138}})));
      equation
        connect(Retha.q_out, inertia.q_out) annotation (
          Line(points = {{112, 17}, {130, 17}}, thickness = 1));
        connect(RightAtrium.q_in[1], TricuspidValve.q_in) annotation (
          Line(points={{-228.14,57.09},{-214,57.09},{-214,58},{-202,58}},          color = {127, 0, 0}, thickness = 0.5));
        connect(TricuspidValve.q_out, RightVentricle.q_in[1]) annotation (
          Line(points={{-176,58},{-166,58},{-166,56.025},{-155.15,56.025}},        color = {127, 0, 0}, thickness = 0.5));
        connect(RightVentricle.q_in[2], PulmonaryValve.q_in) annotation (
          Line(points={{-155.15,57.975},{-143.725,57.975},{-143.725,57},{-132,
              57}},                                                                        color = {127, 0, 0}, thickness = 0.5));
        connect(Epa.q_in[1], Rpp.q_in) annotation (
          Line(points={{-80.14,97.09},{-68,97.09},{-68,98},{-56,98}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rpp.q_out, Epv.q_in[1]) annotation (
          Line(points={{-22,98},{-8,98},{-8,97.09},{2.83,97.09}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Epv.q_in[2], Rlain.q_in) annotation (
          Line(points={{2.83,98.91},{16,98.91},{16,98},{26,98}},          color = {127, 0, 0}, thickness = 0.5));
        connect(LeftAtrium.q_in[1], MitralValve.q_in) annotation (
          Line(points={{87.86,63.09},{102,63.09},{102,64},{114,64}},          color = {127, 0, 0}, thickness = 0.5));
        connect(MitralValve.q_out, LeftVentricle.q_in[1]) annotation (
          Line(points={{140,64},{154,64},{154,66},{163.86,66},{163.86,63.09}},            color = {127, 0, 0}, thickness = 0.5));
        connect(LeftVentricle.q_in[2], AorticValve.q_in) annotation (
          Line(points={{163.86,64.91},{172,64.91},{172,64},{184,64}},          color = {127, 0, 0}, thickness = 0.5));
        connect(inertia.q_in, Eitha.q_in[1]) annotation (
          Line(points={{152,17},{164,17},{164,15.9275},{178.89,15.9275}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Retha.q_in, Eetha.q_in[1]) annotation (
          Line(points={{90,17},{80,17},{80,15.8733},{68.87,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rsart.q_in, Eetha.q_in[2]) annotation (
          Line(points = {{38, 19}, {52, 19}, {52, 17}, {68.87, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Est.q_in[1], Rsart.q_out) annotation (
          Line(points={{-16.12,16.0467},{-3.18,16.0467},{-3.18,19},{10,19}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Est.q_in[2], Rsven.q_in) annotation (
          Line(points = {{-16.12, 17}, {-31.18, 17}, {-31.18, 17}, {-46, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Rsven.q_out, Eethv.q_in[1]) annotation (
          Line(points={{-74,17},{-86,17},{-86,15.8733},{-95.13,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rethv.q_in, Eethv.q_in[2]) annotation (
          Line(points = {{-120, 17}, {-107, 17}, {-107, 17}, {-95.13, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Rethv.q_out, Eithv.q_in[1]) annotation (
          Line(points={{-146,17},{-164,17},{-164,15.8733},{-180.14,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rrain.q_in, Eithv.q_in[2]) annotation (
          Line(points = {{-208, 17}, {-194, 17}, {-194, 17}, {-180.14, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(EithaPressure.port, Eitha.q_in[2]) annotation (
          Line(points={{200,-36},{178.89,-36},{178.89,16.6425}},       color = {127, 0, 0}, thickness = 0.5));
        connect(EethaPressure.port, Eetha.q_in[3]) annotation (
          Line(points={{88,-56},{88,-54},{68.87,-54},{68.87,18.1267}},          color = {127, 0, 0}, thickness = 0.5));
        connect(EstPressure.port, Est.q_in[3]) annotation (
          Line(points={{-2,-52},{-16.12,-52},{-16.12,17.9533}},        color = {127, 0, 0}, thickness = 0.5));
        connect(EethvPressure.port, Eethv.q_in[3]) annotation (
          Line(points={{-82,-52},{-95.13,-52},{-95.13,18.1267}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Eithv.q_in[3],EithvPressure.port)  annotation (
          Line(points={{-180.14,18.1267},{-180.14,-54},{-156,-54}},        color = {127, 0, 0}, thickness = 0.5));
        connect(RAtrialElastance.Et, RightAtrium.elastance) annotation (
          Line(points = {{-202.39, 101.84}, {-202.39, 85.92}, {-225.2, 85.92}, {-225.2, 70.6}}, color = {0, 0, 127}));
        connect(RVentricularElastance.Et, RightVentricle.elastance) annotation (
          Line(points = {{-147.15, 108.91}, {-134, 108.91}, {-134, 76}, {-152, 76}, {-152, 70.5}}, color = {0, 0, 127}));
        connect(LAtrialElastance.Et, LeftAtrium.elastance) annotation (
          Line(points = {{121.61, 107.84}, {121.61, 91.92}, {90.8, 91.92}, {90.8, 76.6}}, color = {0, 0, 127}));
        connect(LVentricularElastance.Et, LeftVentricle.elastance) annotation (
          Line(points = {{203.42, 107.68}, {222, 107.68}, {222, 76.6}, {166.8, 76.6}}, color = {0, 0, 127}));
        connect(PulmonaryValve.q_out, power.q_in) annotation (
          Line(points = {{-106, 57}, {-96, 57}, {-96, 62}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power.q_out, Epa.q_in[2]) annotation (
          Line(points={{-96,82},{-96,98.91},{-80.14,98.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Rrain.q_out, power1.q_in) annotation (
          Line(points = {{-236, 17}, {-254, 17}, {-254, 26}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power1.q_out, RightAtrium.q_in[2]) annotation (
          Line(points={{-254,46},{-254,58.91},{-228.14,58.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(power.power, feedback.u1) annotation (
          Line(points = {{-108, 72}, {-276, 72}, {-276, -56}, {-262, -56}}, color = {0, 0, 127}));
        connect(power1.power, feedback.u2) annotation (
          Line(points = {{-266, 36}, {-268, 36}, {-268, -30}, {-254, -30}, {-254, -48}}, color = {0, 0, 127}));
        connect(feedback.y, rightHeartPower.u) annotation (
          Line(points = {{-245, -56}, {-216, -56}}, color = {0, 0, 127}));
        connect(Rlain.q_out, power2.q_in) annotation (
          Line(points = {{56, 98}, {70, 98}, {70, 92}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power2.q_out, LeftAtrium.q_in[2]) annotation (
          Line(points={{70,72},{70,64.91},{87.86,64.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(AorticValve.q_out, power3.q_in) annotation (
          Line(points = {{208, 64}, {224, 64}, {224, 48}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power3.q_out, Eitha.q_in[3]) annotation (
          Line(points={{224,28},{224,17.3575},{178.89,17.3575}},      color = {127, 0, 0}, thickness = 0.5));
        connect(feedback1.y, leftHeartPower.u) annotation (
          Line(points = {{155, -100}, {184, -100}}, color = {0, 0, 127}));
        connect(power3.power, feedback1.u1) annotation (
          Line(points = {{236, 38}, {244, 38}, {244, -84}, {130, -84}, {130, -100}, {138, -100}}, color = {0, 0, 127}));
        connect(power2.power, feedback1.u2) annotation (
          Line(points = {{82, 82}, {106, 82}, {106, -8}, {118, -8}, {118, -124}, {146, -124}, {146, -108}}, color = {0, 0, 127}));
        connect(Eitha.q_in[4],arterialPressure.port)  annotation (
          Line(points={{178.89,18.0725},{178.89,-66},{208,-66}},       color = {127, 0, 0}, thickness = 0.5));
        connect(uniformNoise.y, randomizedHeartRate.u1)
          annotation (Line(points={{-253,174},{-238,174}}, color={0,0,127}));
        connect(randomizedHeartRate.y, currentHeartReat.u1) annotation (Line(
              points={{-215,168},{-184,168},{-184,156}}, color={0,0,127}));
        connect(currentHeartReat.y, RVentricularElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,134},{-165,134},{-165,118.6}},
              color={0,0,127}));
        connect(currentHeartReat.y, LAtrialElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,136},{99,136},{99,120.8}},
              color={0,0,127}));
        connect(currentHeartReat.y, RAtrialElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,130},{-225,130},{-225,114.8}},
              color={0,0,127}));
        connect(currentHeartReat.y, LVentricularElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,136},{182,136},{182,116.8}},
              color={0,0,127}));
        connect(randomizeHR.y, currentHeartReat.u2) annotation (Line(points={{
                -193,152},{-190,152},{-190,148},{-184,148}}, color={255,0,255}));
        connect(Ecg.HR, currentHeartReat.y) annotation (Line(points={{-90,164},
                {-90,172},{-152,172},{-152,148},{-161,148}}, color={0,0,127}));
        connect(convertSaturationToHeartRate.y[1], randomizedHeartRate.u2)
          annotation (Line(points={{-255,130},{-244,130},{-244,154},{-246,154},
                {-246,162},{-238,162}}, color={0,0,127}));
        connect(saturationInTime.y, convertSaturationToHeartRate.u) annotation (
           Line(points={{-291,128},{-284,128},{-284,130},{-278,130}}, color={0,
                0,127}));
        connect(convertSaturationToHeartRate.y[1], currentHeartReat.u3)
          annotation (Line(points={{-255,130},{-244,130},{-244,134},{-192,134},
                {-192,132},{-184,132},{-184,140}}, color={0,0,127}));
        annotation (
          Diagram(coordinateSystem(extent = {{-280, -140}, {280, 180}}, preserveAspectRatio = false)),
          Icon(coordinateSystem(extent = {{-280, -140}, {280, 180}}, preserveAspectRatio = false), graphics),
          Documentation(info = "<html>
        <p>Model of cardiovascular system using to demonstrate elastic and resistance features of veins and arteries in pulmonary and systemic circulation and influence of cardiac output on it.</p>
        <ul>
        <li>J. A. Goodwin, W. L. van Meurs, C. D. Sa Couto, J. E. W.Beneken, S. A. Graves, A model for educational simulation of infant cardiovascular physiology., Anesthesia and analgesia 99 (6)(2004) 1655&ndash;1664. doi:10.1213/01.ANE.0000134797.52793.AF.</li>
        <li>C. D. Sa Couto, W. L. van Meurs, J. A. Goodwin, P. Andriessen,A Model for Educational Simulation of Neonatal Cardiovascular Pathophysiology, Simulation in Healthcare 1 (Inaugural) (2006) 4&ndash;12.</li>
        <li>W. van Meurs, Modeling and Simulation in Biomedical Engineering: Applications in Cardiorespiratory Physiology, McGraw-Hill Professional, 2011.</li>
        </ul>
        </html>",
                revisions = "<html>
        <ul>
        <li><i>Jul 2015 </i>by Tomas Kulhanek: Created. </li>
        </ul>
        </html>"),
          experiment(StopTime = 10, __Dymola_Algorithm = "Dassl"));
      end testHemodynamics;

      model HemodynamicsRegulatedHR
        extends Physiolibrary.Icons.CardioVascular;
        Physiolibrary.Fluid.Components.ElasticVesselElastance Epa(
        volume_start=0.000106,
        ZeroPressureVolume=5e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,
        Elastance=31064116.267695,                                                                                                                                                            nPorts = 2) annotation (
          Placement(transformation(extent = {{-94, 84}, {-66, 112}})));
        Physiolibrary.Fluid.Components.Resistor Rpp(Resistance(displayUnit = "(mmHg.s)/ml") = 14665462.61565) annotation (
          Placement(transformation(extent = {{-56, 85}, {-22, 111}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Epv(
        volume_start=0.000518,
        ZeroPressureVolume=0.00035,
        ExternalPressure=-533.28954966,
        Elastance=6066168.6273825,                                                                                                                                                              nPorts = 2) annotation (
          Placement(transformation(extent = {{-14, 84}, {20, 112}})));
        Physiolibrary.Fluid.Components.Resistor Rlain(Resistance(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{26, 86}, {56, 110}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance LeftAtrium(
        isExternalPressureAbsolute=false,                                useElastanceInput = true,
        volume_start=9.31e-05,
        ZeroPressureVolume=3e-05,
        ExternalPressure=-533.28954966,                                                                                                                                                           nPorts = 2) annotation (
          Placement(transformation(extent = {{74, 50}, {102, 78}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance LeftVentricle(
        isExternalPressureAbsolute=false,                                   useElastanceInput = true,
        volume_start=0.000144,
        ZeroPressureVolume=6e-05,
        ExternalPressure=-533.28954966,                                                                                                                                                              nPorts = 2) annotation (
          Placement(transformation(extent = {{150, 50}, {178, 78}})));
        Physiolibrary.Fluid.Components.IdealValveResistance AorticValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 1066579.09932) annotation (
          Placement(transformation(extent = {{184, 76}, {208, 52}})));
        Parts.AtrialElastance LAtrialElastance(Tav(displayUnit = "s"), EMIN = 15998686.4898, EMAX = 37330268.4762) annotation (
          Placement(transformation(extent = {{80, 92}, {118, 124}})));
        Parts.VentricularElastance LVentricularElastance(EMIN = 11999014.86735, EMAX = 533289549.66) annotation (
          Placement(transformation(extent = {{164, 88}, {200, 120}})));
        Physiolibrary.Fluid.Components.IdealValveResistance MitralValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(origin = {127, 64}, extent = {{-13, 12}, {13, -12}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eitha(
        ExternalPressure=-533.28954966,                                                                       nPorts = 4,
        volume_start=0.0002,
        ZeroPressureVolume=0.00014,
        Elastance=190651014.00345)                                                                                                                                                                          annotation (
          Placement(transformation(extent = {{168, 6}, {190, 28}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eetha(
        volume_start=0.000526,
        ZeroPressureVolume=0.00037,
        Elastance=74127247.40274,                                                                                                                      nPorts = 3) annotation (
          Placement(transformation(extent = {{56, 4}, {82, 30}})));
        Physiolibrary.Fluid.Components.Inertia inertia(I(displayUnit=
              "mmHg.s2/g") = 226.6480586055, massFlow_start(displayUnit=
              "g/min") = 0.021666666666667)                                                                                                                       annotation (
          Placement(transformation(extent = {{-11, -11}, {11, 11}}, rotation = 180, origin = {141, 17})));
        Physiolibrary.Fluid.Components.Resistor Retha(Resistance(displayUnit = "(mmHg.s)/ml") = 7999343.2449) annotation (
          Placement(transformation(extent = {{90, 6}, {112, 28}})));
        Physiolibrary.Fluid.Components.Resistor Rsart(Resistance(displayUnit = "(mmHg.s)/ml") = 106657909.932) annotation (
          Placement(transformation(extent = {{14, -13}, {-14, 13}}, origin = {24, 19})));
        Physiolibrary.Fluid.Components.Resistor Rsven(Resistance(displayUnit = "(mmHg.s)/ml") = 26664477.483) annotation (
          Placement(transformation(extent = {{14, -13}, {-14, 13}}, origin = {-60, 17})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Est(volume_start = 0.283e-3, ZeroPressureVolume = 0.185e-3, Elastance = 34930465.50273, nPorts = 3) annotation (
          Placement(transformation(extent = {{-28, 6}, {-4, 28}})));
        Physiolibrary.Fluid.Components.Resistor Rethv(Resistance(displayUnit = "(mmHg.s)/ml") = 11999014.86735) annotation (
          Placement(transformation(extent = {{-120, 4}, {-146, 30}})));
        Physiolibrary.Fluid.Components.Resistor Rrain(Resistance(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{-208, 4}, {-236, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eithv(
        volume_start=0.00148,
        ZeroPressureVolume=0.00119,
        ExternalPressure=-533.28954966,
        Elastance=2426467.450953,                                                                                                                                                               nPorts = 3) annotation (
          Placement(transformation(extent = {{-194, 4}, {-166, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance Eethv(
        volume_start=0.00153,
        ZeroPressureVolume=0.001,
        Elastance=2253148.3473135,                                                                                                                  nPorts = 3) annotation (
          Placement(transformation(extent = {{-108, 4}, {-82, 30}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance RightAtrium(
        isExternalPressureAbsolute=false,                                 useElastanceInput = true,
        volume_start=0.000135,
        ZeroPressureVolume=3e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,                                                                                                                                      nPorts = 2) annotation (
          Placement(transformation(extent = {{-242, 44}, {-214, 72}})));
        Physiolibrary.Fluid.Components.ElasticVesselElastance RightVentricle(
        isExternalPressureAbsolute=false,                                    useElastanceInput = true,
        volume_start=0.000131,
        ZeroPressureVolume=4e-05,
        ExternalPressure(displayUnit="mmHg") = -533.28954966,                                                                                                                                         nPorts = 2) annotation (
          Placement(transformation(extent = {{-170, 42}, {-140, 72}})));
        Physiolibrary.Fluid.Components.IdealValveResistance PulmonaryValve(_Roff(displayUnit = "g/(mmHg.s)") = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(extent = {{-132, 70}, {-106, 44}})));
        Parts.AtrialElastance RAtrialElastance(EMIN = 6666119.37075, EMAX = 19998358.11225) annotation (
          Placement(transformation(extent = {{-244, 86}, {-206, 118}})));
        Parts.VentricularElastance RVentricularElastance(EMIN = 7599376.082655, EMAX = 65327969.83335) annotation (
          Placement(transformation(extent = {{-180, 88}, {-150, 122}})));
        Physiolibrary.Fluid.Components.IdealValveResistance TricuspidValve(_Roff = Modelica.Constants.inf, _Ron(displayUnit = "(mmHg.s)/ml") = 399967.162245) annotation (
          Placement(transformation(origin = {-189, 58}, extent = {{-13, 12}, {13, -12}})));
        inner Modelica.Fluid.System system(p_ambient(displayUnit = "mmHg") = 101325.0144354, T_ambient = 310.15) annotation (
          Placement(transformation(extent = {{-48, 144}, {-28, 164}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EithaPressure
          annotation (Placement(transformation(extent={{190,-36},{210,-16}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EethaPressure
          annotation (Placement(transformation(extent={{78,-56},{98,-36}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EstPressure
          annotation (Placement(transformation(extent={{-12,-52},{8,-32}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EethvPressure
          annotation (Placement(transformation(extent={{-92,-52},{-72,-32}})));
        Physiolibrary.Fluid.Sensors.PressureMeasure EithvPressure annotation (
            Placement(transformation(extent={{-166,-54},{-146,-34}})));
        Physiolibrary.Fluid.Sensors.Power power annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-96,72})));
        Physiolibrary.Fluid.Sensors.Power power1 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-254,36})));
        Modelica.Blocks.Math.Mean rightHeartPower(f(displayUnit = "1/min") = 1.2) annotation (
          Placement(transformation(extent = {{-214, -66}, {-194, -46}})));
        Modelica.Blocks.Math.Feedback feedback annotation (
          Placement(transformation(extent = {{-264, -46}, {-244, -66}})));
        Physiolibrary.Fluid.Sensors.Power power2 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,82})));
        Physiolibrary.Fluid.Sensors.Power power3 annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={224,38})));
        Modelica.Blocks.Math.Mean leftHeartPower(f(displayUnit = "1/min") = 1.2) annotation (
          Placement(transformation(extent = {{186, -110}, {206, -90}})));
        Modelica.Blocks.Math.Feedback feedback1 annotation (
          Placement(transformation(extent = {{136, -110}, {156, -90}})));
        Physiolibrary.Fluid.Sensors.Sphygmomanometer arterialPressure(
            MeasurementTime(displayUnit="s") = 60/72)
          annotation (Placement(transformation(extent={{198,-66},{218,-46}})));
        Modelica.Blocks.Noise.UniformNoise uniformNoise(
          samplePeriod=6,
          y_min=-0.03,
          y_max=0.03) annotation (Placement(transformation(extent={{-274,164},{
                  -254,184}})));
        Modelica.Blocks.Math.Add randomizedHeartRate annotation (Placement(
              transformation(extent={{-236,158},{-216,178}})));
        Modelica.Blocks.Logical.Switch currentHeartReat annotation (Placement(
              transformation(extent={{-182,138},{-162,158}})));
        Modelica.Blocks.Sources.BooleanConstant randomizeHR annotation (
            Placement(transformation(extent={{-214,142},{-194,162}})));
        Parts.ECGGenerator Ecg
          annotation (Placement(transformation(extent={{-100,146},{-80,166}})));
        Physiolibrary.Types.Constants.FractionConst sO2(k=0.98)   annotation (
            Placement(transformation(extent={{-340,128},{-332,136}})));
        Modelica.Blocks.Tables.CombiTable1Ds convertSaturationToHeartRate(table
            =[0,0.03; 0.55,0.03; 0.6,180/60; 0.7,120/60; 0.8,95/60; 0.9,72/60],
            extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
          annotation (Placement(transformation(extent={{-310,120},{-290,140}})));
        Modelica.Blocks.Math.Add add annotation (Placement(transformation(
                extent={{-268,130},{-248,150}})));
        Modelica.Blocks.Sources.Constant HRAdd(k=0) annotation (Placement(
              transformation(extent={{-320,148},{-300,168}})));
      equation
        connect(Retha.q_out, inertia.q_out) annotation (
          Line(points = {{112, 17}, {130, 17}}, thickness = 1));
        connect(RightAtrium.q_in[1], TricuspidValve.q_in) annotation (
          Line(points={{-228.14,57.09},{-214,57.09},{-214,58},{-202,58}},          color = {127, 0, 0}, thickness = 0.5));
        connect(TricuspidValve.q_out, RightVentricle.q_in[1]) annotation (
          Line(points={{-176,58},{-166,58},{-166,56.025},{-155.15,56.025}},        color = {127, 0, 0}, thickness = 0.5));
        connect(RightVentricle.q_in[2], PulmonaryValve.q_in) annotation (
          Line(points={{-155.15,57.975},{-143.725,57.975},{-143.725,57},{-132,
              57}},                                                                        color = {127, 0, 0}, thickness = 0.5));
        connect(Epa.q_in[1], Rpp.q_in) annotation (
          Line(points={{-80.14,97.09},{-68,97.09},{-68,98},{-56,98}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rpp.q_out, Epv.q_in[1]) annotation (
          Line(points={{-22,98},{-8,98},{-8,97.09},{2.83,97.09}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Epv.q_in[2], Rlain.q_in) annotation (
          Line(points={{2.83,98.91},{16,98.91},{16,98},{26,98}},          color = {127, 0, 0}, thickness = 0.5));
        connect(LeftAtrium.q_in[1], MitralValve.q_in) annotation (
          Line(points={{87.86,63.09},{102,63.09},{102,64},{114,64}},          color = {127, 0, 0}, thickness = 0.5));
        connect(MitralValve.q_out, LeftVentricle.q_in[1]) annotation (
          Line(points={{140,64},{154,64},{154,66},{163.86,66},{163.86,63.09}},            color = {127, 0, 0}, thickness = 0.5));
        connect(LeftVentricle.q_in[2], AorticValve.q_in) annotation (
          Line(points={{163.86,64.91},{172,64.91},{172,64},{184,64}},          color = {127, 0, 0}, thickness = 0.5));
        connect(inertia.q_in, Eitha.q_in[1]) annotation (
          Line(points={{152,17},{164,17},{164,15.9275},{178.89,15.9275}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Retha.q_in, Eetha.q_in[1]) annotation (
          Line(points={{90,17},{80,17},{80,15.8733},{68.87,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rsart.q_in, Eetha.q_in[2]) annotation (
          Line(points = {{38, 19}, {52, 19}, {52, 17}, {68.87, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Est.q_in[1], Rsart.q_out) annotation (
          Line(points={{-16.12,16.0467},{-3.18,16.0467},{-3.18,19},{10,19}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Est.q_in[2], Rsven.q_in) annotation (
          Line(points = {{-16.12, 17}, {-31.18, 17}, {-31.18, 17}, {-46, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Rsven.q_out, Eethv.q_in[1]) annotation (
          Line(points={{-74,17},{-86,17},{-86,15.8733},{-95.13,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rethv.q_in, Eethv.q_in[2]) annotation (
          Line(points = {{-120, 17}, {-107, 17}, {-107, 17}, {-95.13, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(Rethv.q_out, Eithv.q_in[1]) annotation (
          Line(points={{-146,17},{-164,17},{-164,15.8733},{-180.14,15.8733}},          color = {127, 0, 0}, thickness = 0.5));
        connect(Rrain.q_in, Eithv.q_in[2]) annotation (
          Line(points = {{-208, 17}, {-194, 17}, {-194, 17}, {-180.14, 17}}, color = {127, 0, 0}, thickness = 0.5));
        connect(EithaPressure.port, Eitha.q_in[2]) annotation (
          Line(points={{200,-36},{178.89,-36},{178.89,16.6425}},       color = {127, 0, 0}, thickness = 0.5));
        connect(EethaPressure.port, Eetha.q_in[3]) annotation (
          Line(points={{88,-56},{88,-54},{68.87,-54},{68.87,18.1267}},          color = {127, 0, 0}, thickness = 0.5));
        connect(EstPressure.port, Est.q_in[3]) annotation (
          Line(points={{-2,-52},{-16.12,-52},{-16.12,17.9533}},        color = {127, 0, 0}, thickness = 0.5));
        connect(EethvPressure.port, Eethv.q_in[3]) annotation (
          Line(points={{-82,-52},{-95.13,-52},{-95.13,18.1267}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Eithv.q_in[3],EithvPressure.port)  annotation (
          Line(points={{-180.14,18.1267},{-180.14,-54},{-156,-54}},        color = {127, 0, 0}, thickness = 0.5));
        connect(RAtrialElastance.Et, RightAtrium.elastance) annotation (
          Line(points = {{-202.39, 101.84}, {-202.39, 85.92}, {-225.2, 85.92}, {-225.2, 70.6}}, color = {0, 0, 127}));
        connect(RVentricularElastance.Et, RightVentricle.elastance) annotation (
          Line(points = {{-147.15, 108.91}, {-134, 108.91}, {-134, 76}, {-152, 76}, {-152, 70.5}}, color = {0, 0, 127}));
        connect(LAtrialElastance.Et, LeftAtrium.elastance) annotation (
          Line(points = {{121.61, 107.84}, {121.61, 91.92}, {90.8, 91.92}, {90.8, 76.6}}, color = {0, 0, 127}));
        connect(LVentricularElastance.Et, LeftVentricle.elastance) annotation (
          Line(points = {{203.42, 107.68}, {222, 107.68}, {222, 76.6}, {166.8, 76.6}}, color = {0, 0, 127}));
        connect(PulmonaryValve.q_out, power.q_in) annotation (
          Line(points = {{-106, 57}, {-96, 57}, {-96, 62}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power.q_out, Epa.q_in[2]) annotation (
          Line(points={{-96,82},{-96,98.91},{-80.14,98.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(Rrain.q_out, power1.q_in) annotation (
          Line(points = {{-236, 17}, {-254, 17}, {-254, 26}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power1.q_out, RightAtrium.q_in[2]) annotation (
          Line(points={{-254,46},{-254,58.91},{-228.14,58.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(power.power, feedback.u1) annotation (
          Line(points = {{-108, 72}, {-276, 72}, {-276, -56}, {-262, -56}}, color = {0, 0, 127}));
        connect(power1.power, feedback.u2) annotation (
          Line(points = {{-266, 36}, {-268, 36}, {-268, -30}, {-254, -30}, {-254, -48}}, color = {0, 0, 127}));
        connect(feedback.y, rightHeartPower.u) annotation (
          Line(points = {{-245, -56}, {-216, -56}}, color = {0, 0, 127}));
        connect(Rlain.q_out, power2.q_in) annotation (
          Line(points = {{56, 98}, {70, 98}, {70, 92}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power2.q_out, LeftAtrium.q_in[2]) annotation (
          Line(points={{70,72},{70,64.91},{87.86,64.91}},        color = {127, 0, 0}, thickness = 0.5));
        connect(AorticValve.q_out, power3.q_in) annotation (
          Line(points = {{208, 64}, {224, 64}, {224, 48}}, color = {127, 0, 0}, thickness = 0.5));
        connect(power3.q_out, Eitha.q_in[3]) annotation (
          Line(points={{224,28},{224,17.3575},{178.89,17.3575}},      color = {127, 0, 0}, thickness = 0.5));
        connect(feedback1.y, leftHeartPower.u) annotation (
          Line(points = {{155, -100}, {184, -100}}, color = {0, 0, 127}));
        connect(power3.power, feedback1.u1) annotation (
          Line(points = {{236, 38}, {244, 38}, {244, -84}, {130, -84}, {130, -100}, {138, -100}}, color = {0, 0, 127}));
        connect(power2.power, feedback1.u2) annotation (
          Line(points = {{82, 82}, {106, 82}, {106, -8}, {118, -8}, {118, -124}, {146, -124}, {146, -108}}, color = {0, 0, 127}));
        connect(Eitha.q_in[4],arterialPressure.port)  annotation (
          Line(points={{178.89,18.0725},{178.89,-66},{208,-66}},       color = {127, 0, 0}, thickness = 0.5));
        connect(uniformNoise.y, randomizedHeartRate.u1)
          annotation (Line(points={{-253,174},{-238,174}}, color={0,0,127}));
        connect(randomizedHeartRate.y, currentHeartReat.u1) annotation (Line(
              points={{-215,168},{-184,168},{-184,156}}, color={0,0,127}));
        connect(currentHeartReat.y, RVentricularElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,134},{-165,134},{-165,118.6}},
              color={0,0,127}));
        connect(currentHeartReat.y, LAtrialElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,136},{99,136},{99,120.8}},
              color={0,0,127}));
        connect(currentHeartReat.y, RAtrialElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,130},{-225,130},{-225,114.8}},
              color={0,0,127}));
        connect(currentHeartReat.y, LVentricularElastance.HR) annotation (Line(
              points={{-161,148},{-152,148},{-152,136},{182,136},{182,116.8}},
              color={0,0,127}));
        connect(randomizeHR.y, currentHeartReat.u2) annotation (Line(points={{
                -193,152},{-190,152},{-190,148},{-184,148}}, color={255,0,255}));
        connect(Ecg.HR, currentHeartReat.y) annotation (Line(points={{-90,164},
                {-90,172},{-152,172},{-152,148},{-161,148}}, color={0,0,127}));
        connect(sO2.y, convertSaturationToHeartRate.u)
          annotation (Line(points={{-331,132},{-320,132},{-320,130},{-312,130}},
                                                           color={0,0,127}));
        connect(convertSaturationToHeartRate.y[1], add.u2) annotation (Line(
              points={{-289,130},{-278,130},{-278,134},{-270,134}}, color={0,0,
                127}));
        connect(add.y, randomizedHeartRate.u2) annotation (Line(points={{-247,
                140},{-242,140},{-242,124},{-282,124},{-282,158},{-248,158},{
                -248,162},{-238,162}}, color={0,0,127}));
        connect(add.y, currentHeartReat.u3) annotation (Line(points={{-247,140},
                {-242,140},{-242,134},{-192,134},{-192,132},{-184,132},{-184,
                140}}, color={0,0,127}));
        connect(HRAdd.y, add.u1) annotation (Line(points={{-299,158},{-284,158},
                {-284,146},{-270,146}}, color={0,0,127}));
        annotation (
          Diagram(coordinateSystem(extent = {{-280, -140}, {280, 180}}, preserveAspectRatio = false)),
          Icon(coordinateSystem(extent = {{-280, -140}, {280, 180}}, preserveAspectRatio = false), graphics),
          Documentation(info = "<html>
        <p>Model of cardiovascular system using to demonstrate elastic and resistance features of veins and arteries in pulmonary and systemic circulation and influence of cardiac output on it.</p>
        <ul>
        <li>J. A. Goodwin, W. L. van Meurs, C. D. Sa Couto, J. E. W.Beneken, S. A. Graves, A model for educational simulation of infant cardiovascular physiology., Anesthesia and analgesia 99 (6)(2004) 1655&ndash;1664. doi:10.1213/01.ANE.0000134797.52793.AF.</li>
        <li>C. D. Sa Couto, W. L. van Meurs, J. A. Goodwin, P. Andriessen,A Model for Educational Simulation of Neonatal Cardiovascular Pathophysiology, Simulation in Healthcare 1 (Inaugural) (2006) 4&ndash;12.</li>
        <li>W. van Meurs, Modeling and Simulation in Biomedical Engineering: Applications in Cardiorespiratory Physiology, McGraw-Hill Professional, 2011.</li>
        </ul>
        </html>",
                revisions = "<html>
        <ul>
        <li><i>Jul 2015 </i>by Tomas Kulhanek: Created. </li>
        </ul>
        </html>"),
          experiment(StopTime = 10, __Dymola_Algorithm = "Dassl"));
      end HemodynamicsRegulatedHR;
      annotation (
        Documentation(info = "<html>
	</html>"));
    end MeursModel2011;

    model LungVentilatorSCMV2
      extends Physiolibrary.Icons.RespiratoryCenter;
      replaceable package Air = Physiolibrary.Media.Air;
      parameter Physiolibrary.Types.HydraulicResistance TotalResistance = 147099.75 "Total lungs pathways conductance";
      parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(displayUnit = "ml/mmHg") = 6.0004926067653e-07 "Total lungs compliance";
      parameter Physiolibrary.Types.Volume TV=0.0005                                               "Tidal volume [m3]";
      parameter Modelica.Units.SI.Volume LungsAirVolume_initial=FunctionalResidualCapacity;
      parameter Modelica.Units.SI.Volume FunctionalResidualCapacity=0.00231 "Functional residual capacity";
      parameter Physiolibrary.Types.Frequency RR=0.28333333333333 "Respiration rate [s-1]";

      Modelica.Blocks.Math.Add add1(k2 = +1) annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin={-22,30})));
      Physiolibrary.Types.Constants.VolumeConst deadspace_c(k(displayUnit = "l") = 0.00035) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin={-4,-6})));
      Modelica.Blocks.Math.Division respirationRate annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin={-24,74})));
      Modelica.Blocks.Math.Product RR_multiply annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 180, origin={-35,95})));
      Physiolibrary.Types.Constants.VolumeConst base_deadspace_volume(k = DV) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 180, origin={-4,96})));
      Physiolibrary.Types.Constants.VolumeConst volume(k(displayUnit = "l") = 0.0023) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin={-6,30})));
      Modelica.Blocks.Math.Min tidalVolume annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 90, origin={-19,53})));
      Modelica.Blocks.Math.Product product4 annotation (
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 90, origin={-38,-6})));
      Physiolibrary.Types.Constants.FrequencyConst m(k = 0.505) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 180, origin={-12,-32})));
      parameter Physiolibrary.Types.Pressure pCO2_ZeroVentilation(displayUnit = "kPa") = 4800 "Long term adaptation for pCO2 during acidosis/alcalosis";
      parameter Physiolibrary.Types.Volume DV = 0.00015 "Deadspace volume";
      Modelica.Blocks.Interfaces.RealOutput deadSpaceVentilation annotation (
        Placement(transformation(rotation = 0, extent={{12,45.5},{32,70.5}}),       iconTransformation(extent = {{-10, -12.5}, {10, 12.5}}, rotation = 270, origin = {60, -108})));
      Modelica.Blocks.Interfaces.RealOutput ventilation annotation (
        Placement(transformation(rotation = 0, extent={{-106,116.5},{-86,141.5}}),   iconTransformation(extent = {{-10, -12.5}, {10, 12.5}}, rotation = 270, origin = {4, -107})));
      Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
          Medium = Air)                                                                      "External environment" annotation (
        Placement(transformation(extent={{-74,-106},{-54,-86}})));
      Physiolibrary.Fluid.Sensors.FlowMeasure
                          flowMeasure(redeclare package Medium = Air) annotation (
        Placement(transformation(extent={{48,-112},{68,-92}})));
      Physiolibrary.Fluid.Components.Resistor
                          resistor(redeclare package Medium = Air, Resistance=
            TotalResistance)                                                                      annotation (
        Placement(transformation(extent={{78,-112},{98,-92}})));
      Physiolibrary.Types.Constants.FrequencyConst
                                     frequency(k=RR)                annotation (
        Placement(transformation(extent={{-90,-60},{-74,-46}})));
      Physiolibrary.Fluid.Components.VolumePump volumePump(redeclare package
          Medium =
            Air, useSolutionFlowInput=true)
        annotation (Placement(transformation(extent={{-7,-112},{13,-92}})));
      VolumeController ventilatorSCMV(peakvolume=TV, maxflowrate(displayUnit="l/s"))
        annotation (Placement(transformation(extent={{-1,-84},{19,-64}})));
      Physiolibrary.Fluid.Components.Conductor expiration(redeclare package
          Medium =
            Air, useConductanceInput=true) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={47,-58})));
      Physiolibrary.Fluid.Components.IdealValve idealValve(redeclare package
          Medium =
            Air) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={69,12})));
      Physiolibrary.Fluid.Sources.PressureSource environment1(redeclare package
          Medium = Air)                                                                      "External environment" annotation (
        Placement(transformation(extent={{10,2},{30,22}})));
      inner Modelica.Fluid.System system(T_ambient=310.15)   annotation (
        Placement(transformation(extent={{40,120},{60,140}})));
      Physiolibrary.Fluid.Components.ElasticVessel lungs(
        redeclare package Medium = Air,
        use_mass_start=false,
        volume_start=LungsAirVolume_initial,
        massFractions_start=Air.reference_X,
        ZeroPressureVolume=FunctionalResidualCapacity,
        Compliance=TotalCompliance,
        useExternalPressureInput=false,
        nPorts=2)                                                                                                                                                            "Lungs" annotation (
        Placement(transformation(extent={{112,-112},{132,-92}})));
      Physiolibrary.Fluid.Sensors.PressureMeasure lungsPressureMeasure(
          redeclare package Medium = Air)                                                              "Lungs pressure" annotation (
        Placement(transformation(extent={{148,-104},{168,-84}})));
      Modelica.Blocks.Continuous.Filter filter(
        analogFilter=Modelica.Blocks.Types.AnalogFilter.Butterworth,
        filterType=Modelica.Blocks.Types.FilterType.LowPass,
        order=4,
        f_cut=0.5)
        annotation (Placement(transformation(extent={{-50,-62},{-30,-42}})));
      Physiolibrary.Types.Constants.HydraulicConductanceConst expirationConductance(k(
            displayUnit="m3/(Pa.s)") = 1e-5)
        annotation (Placement(transformation(extent={{11,-24},{19,-16}})));
      Modelica.Blocks.Math.Product product1
        annotation (Placement(transformation(extent={{14,-54},{34,-34}})));
    equation
      connect(deadspace_c.y, add1.u2) annotation (
        Line(points={{-4,-1},{-4,6},{-18.4,6},{-18.4,22.8}},             color = {0, 0, 127}));
      connect(base_deadspace_volume.y, RR_multiply.u2) annotation (
        Line(points={{-9,96},{-9,99.2},{-26.6,99.2}},       color = {0, 0, 127}));
      connect(respirationRate.y, RR_multiply.u1) annotation (
        Line(points={{-24,80.6},{-12,80.6},{-12,88},{-20,88},{-20,90.8},{-26.6,90.8}},        color = {0, 0, 127}));
      connect(tidalVolume.u1, add1.y) annotation (
        Line(points={{-22,47},{-22,36.6}},   color = {0, 0, 127}));
      connect(tidalVolume.y, respirationRate.u2) annotation (
        Line(points={{-19,58.5},{-20.4,58.5},{-20.4,66.8}},     color = {0, 0, 127}));
      connect(volume.y, tidalVolume.u2) annotation (
        Line(points={{-6,35},{-6,47},{-16,47}},      color = {0, 0, 127}));
      connect(product4.y, add1.u1) annotation (
        Line(points={{-38,0.6},{-38,22.8},{-25.6,22.8}},       color = {0, 0, 127}));
      connect(m.y, product4.u2) annotation (
        Line(points={{-17,-32},{-34.4,-32},{-34.4,-13.2}},                           color = {0, 0, 127}));
      connect(deadSpaceVentilation, RR_multiply.y) annotation (
        Line(points={{22,58},{22,108},{-48,108},{-48,95},{-42.7,95}},                  color = {0, 0, 127}));
      connect(flowMeasure.q_out, resistor.q_in) annotation (Line(
          points={{68,-102},{78,-102}},
          color={127,0,0},
          thickness=0.5));
      connect(volumePump.q_out, flowMeasure.q_in) annotation (Line(
          points={{13,-102},{48,-102}},
          color={127,0,0},
          thickness=0.5));
      connect(ventilatorSCMV.volumeflowrate, volumePump.solutionFlow)
        annotation (Line(points={{3,-82.6},{3,-95}}, color={0,0,127}));
      connect(environment1.y, idealValve.q_out) annotation (Line(
          points={{30,12},{59,12}},
          color={127,0,0},
          thickness=0.5));
      connect(idealValve.q_in, expiration.q_out) annotation (Line(
          points={{79,12},{86,12},{86,-42},{47,-42},{47,-48}},
          color={127,0,0},
          thickness=0.5));
      connect(expiration.q_in, flowMeasure.q_in) annotation (Line(
          points={{47,-68},{47,-86},{42,-86},{42,-102},{48,-102}},
          color={127,0,0},
          thickness=0.5));
      connect(product4.u1, ventilation) annotation (Line(points={{-41.6,-13.2},
              {-41.6,-18},{-96,-18},{-96,129}}, color={0,0,127}));
      connect(respirationRate.u1, ventilation) annotation (Line(points={{-27.6,
              66.8},{-27.6,62},{-96,62},{-96,129}}, color={0,0,127}));
      connect(lungs.q_in[1], resistor.q_out) annotation (Line(
          points={{121.9,-102.65},{110,-102.65},{110,-102},{98,-102}},
          color={127,0,0},
          thickness=0.5));
      connect(lungs.fluidVolume, ventilatorSCMV.volume) annotation (Line(points
            ={{132,-110},{140,-110},{140,-74},{16.6,-74}}, color={0,0,127}));
      connect(environment.y, volumePump.q_in) annotation (Line(
          points={{-54,-96},{-14,-96},{-14,-102},{-7,-102}},
          color={127,0,0},
          thickness=0.5));
      connect(ventilatorSCMV.meanventilation, ventilation) annotation (Line(
            points={{-1,-70.4},{-96,-70.4},{-96,129}}, color={0,0,127}));
      connect(frequency.y, filter.u) annotation (Line(points={{-72,-53},{-72,
              -54},{-60,-54},{-60,-52},{-52,-52}}, color={0,0,127}));
      connect(filter.y, ventilatorSCMV.frequency) annotation (Line(points={{-29,
              -52},{2.8,-52},{2.8,-68}}, color={0,0,127}));
      connect(ventilatorSCMV.outflowopen, product1.u2) annotation (Line(points=
              {{-0.7,-74.1},{-12,-74.1},{-12,-50},{12,-50}}, color={0,0,127}));
      connect(product1.y, expiration.cond) annotation (Line(points={{35,-44},{
              35,-28},{41,-28},{41,-58}}, color={0,0,127}));
      connect(product1.u1, expirationConductance.y) annotation (Line(points={{
              12,-38},{4,-38},{4,-12},{24,-12},{24,-20},{20,-20}}, color={0,0,
              127}));

      connect(lungsPressureMeasure.port, lungs.q_in[2]) annotation (Line(points
            ={{158,-104},{158,-138},{121.9,-138},{121.9,-101.35}}, color={0,127,
              255}));
      annotation (
        Diagram(coordinateSystem(extent={{-180,-140},{180,160}})),                                                                   Icon(
            coordinateSystem(extent={{-180,-140},{180,160}})),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
    end LungVentilatorSCMV2;

  end BloodGasesTransport;
  annotation (uses(
      Modelica(version="4.0.0"),
      Chemical(version="1.4.0"),
      Physiolibrary(version="3.0.0-beta1")),   version="1");
end modelECMORespiratoryVR;
