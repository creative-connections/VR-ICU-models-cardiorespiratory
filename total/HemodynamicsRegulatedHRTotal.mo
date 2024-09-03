package Chemical "Physical Chemistry"
 extends Modelica.Icons.Package;

  package Substances "Definitions of substances"
      extends Modelica.Icons.Package;

    record Water_liquid "H2O(l) with self-clustering"
     extends Chemical.Interfaces.Incompressible.SubstanceData(
        MolarWeight=0.018015,
        DfH=-285830,
        DfG=-227230,
        Cp=75.3,
        SelfClustering = true,
        SelfClustering_dH = -81.6348,
        SelfClustering_dS = 32.845554,
          References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
                                      //-77.95928,

    /*  SelfClustering_dH = -81.6348,
    SelfClustering_dS = 32.8344,
*/

      // S=(0 + Modelica.Constants.R*(273.15+25)*log(55.345/0.95-1))/(273.15+25),
      // SelfClustering_dS = (SelfClustering_dH + Modelica.Constants.R*(273.15+25)*log((55.345-1)/1))/(273.15+25),
      annotation (preferredView = "info", Documentation(info="<html>
<p><span style=\"font-family: Courier New;\">Even the tabulated formation Gibbs energy is DfG=-237190 there is another values because of water self-clustering. </span></p>
<p><br><span style=\"font-family: Courier New;\">New reported values for free water molecule in solution is calculated from water dissociation reaction.</span></p>
<p><span style=\"font-family: Courier New;\">&nbsp;&nbsp;&nbsp;&nbsp;</span></p>
</html>"));
    end Water_liquid;

    record Proton_aqueous "H+(aq)"
     extends Chemical.Interfaces.Incompressible.SubstanceData(
        MolarWeight=0.001007,
        z=1,
        DfH=0,
        DfG=0,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Proton_aqueous;

    record Oxygen_aqueous "O2(aq)"
     extends Chemical.Interfaces.Incompressible.SubstanceData(
        MolarWeight=0.032,
        DfH=-11700,
        DfG=16320,
        References={
            "http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf, https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"});
          annotation (preferredView = "info");
    end Oxygen_aqueous;

    record Hydroxide_aqueous "OH-(aq)"
     extends Chemical.Interfaces.Incompressible.SubstanceData(
        MolarWeight=0.017006,
        z=-1,
        DfH=-229940,
        DfG=-157300,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Hydroxide_aqueous;

    record Hydrogen_aqueous "H2(aq)"
     extends Chemical.Interfaces.Incompressible.SubstanceData(
        MolarWeight=0.00201588,
        z=0,
        DfH=-4157,
        DfG=17740,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=10#Solubility"});
      annotation (preferredView = "info");
    end Hydrogen_aqueous;
  end Substances;

  package Interfaces "Chemical interfaces"
    extends Modelica.Icons.InterfacesPackage;

    connector SubstancePort
    "Electro-chemical potential and molar change of the substance in the solution"

    Modelica.Units.SI.ChemicalPotential u
      "Electro-chemical potential of the substance in the solution";

    flow Modelica.Units.SI.MolarFlowRate q
      "Molar change of the substance";

      //with molar flow of substance heat energy is changing also..
    stream Modelica.Units.SI.MolarEnthalpy h_outflow
      "Outgoing molar enthalphy";

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg;(T) .. chemical potential of pure substance</p>
<p>DfG(T) .. free Gibbs energy of formation of the substance at current temperature T. </p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
    end SubstancePort;

    connector SubstancePort_a
    "Electro-chemical potential and molar flow of the substance in the solution"
      extends SubstancePort;

    annotation (
        defaultComponentName="port_a",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={158,66,200}),Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={158,66,200},
            fillColor={158,66,200},
            fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid,
              lineThickness=1),
       Text(extent = {{-160,110},{40,50}}, lineColor={172,72,218},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with internal definition of the substance inside the component. </p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstancePort_a;

    partial package StateOfMatter "Abstract package for all state of matters"


     replaceable partial record SubstanceData
        "Definition data of the chemical substance"

     end SubstanceData;


     replaceable function activityCoefficient
      "Return activity coefficient of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

        output Real activityCoefficient "Activity Coefficient";
     end activityCoefficient;

     replaceable function chargeNumberOfIon
      "Return charge number of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.ChargeNumberOfIon chargeNumberOfIon
        "Charge number of ion";
     end chargeNumberOfIon;

     replaceable function molarEnthalpyElectroneutral
      "Molar enthalpy of the substance in electroneutral solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.MolarEnthalpy molarEnthalpyElectroneutral
        "Molar enthalpy";
     end molarEnthalpyElectroneutral;

     function molarEnthalpy
      "Molar enthalpy of the substance with electric potential dependence"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.MolarEnthalpy molarEnthalpy
        "Molar enthalpy";
     algorithm
        molarEnthalpy := molarEnthalpyElectroneutral(substanceData,T,p,v,I) +
             Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
        annotation (Inline=true, smoothOrder=2);
     end molarEnthalpy;

     replaceable function molarEntropyPure
      "Molar entropy of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.MolarEntropy molarEntropyPure
        "Molar entropy of the pure substance";
     end molarEntropyPure;

      function molarEntropy "Molar entropy of the substance in the solution"
            extends Modelica.Icons.Function;
      input Modelica.Units.SI.ChemicalPotential u
        "Electro-chemical potential of the substance";
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.MolarEntropy molarEntropy "Molar entropy";
      algorithm
          molarEntropy :=  (u - molarEnthalpy(substanceData,T,p,v,I))/T;
      end molarEntropy;

     function chemicalPotentialPure "Chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.ChemicalPotential chemicalPotentialPure
        "Base chemical potential";
     algorithm
         chemicalPotentialPure :=  molarEnthalpyElectroneutral(substanceData,T,p,v,I) - T*molarEntropyPure(substanceData,T,p,v,I);
     end chemicalPotentialPure;

     function electroChemicalPotentialPure
      "Electro-chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.ChemicalPotential
        electroChemicalPotentialPure "Base electro-chemical potential";
     algorithm
      electroChemicalPotentialPure := chemicalPotentialPure(
           substanceData,
           T,
           p,
           v,
           I) + Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
     end electroChemicalPotentialPure;

     replaceable function molarVolumePure "Molar volume of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
     end molarVolumePure;

     function molarVolumeExcess
      "Excess molar volume of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.MolarVolume molarVolumeExcess
        "Excess molar volume of the substance in the solution";
     algorithm
        molarVolumeExcess := molarVolumePure(substanceData,T,p,v,I)*
           log(activityCoefficient(substanceData,T,p,v,I)); //zero if activityCoefficient==1
        annotation (Inline=true, smoothOrder=2);
     end molarVolumeExcess;

     replaceable function molarVolume "Molar volume of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.MolarVolume molarVolume "Molar volume";
     algorithm
      molarVolume :=molarVolumePure(
           substanceData,
           T,
           p,
           v,
           I) + molarVolumeExcess(
           substanceData,
           T,
           p,
           v,
           I);
        annotation (Inline=true, smoothOrder=2);
     end molarVolume;

     replaceable function molarHeatCapacityCp
      "Molar heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.MolarHeatCapacity molarHeatCapacityCp
        "Molar heat capacity at constant pressure";
     end molarHeatCapacityCp;

     replaceable function molarMassOfBaseMolecule
        "Molar mass of base molecule of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarMass molarMass "Molar mass";
     end molarMassOfBaseMolecule;

     replaceable function selfClustering "returns true if substance molecules are joining together to clusters"
         extends Modelica.Icons.Function;
            input SubstanceData substanceData "Data record of substance";
            output Boolean selfClustering;
     algorithm
       selfClustering:=false;
     end selfClustering;



     replaceable function selfClusteringBondEnthalpy
      "Enthalpy of joining two base molecules of the substance together to cluster"
         extends Modelica.Icons.Function;
            input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
     algorithm
       selfClusteringEnthalpy:=0;
     end selfClusteringBondEnthalpy;

     replaceable function selfClusteringBondEntropy
      "Entropy of joining two base molecules of the substance together to cluster"
         extends Modelica.Icons.Function;
            input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
     algorithm
       selfClusteringEntropy:=0;
     end selfClusteringBondEntropy;

     replaceable function selfClusteringBondVolume
         extends Modelica.Icons.Function;
            input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarVolume selfClusteringBondVolume;
     algorithm
       selfClusteringBondVolume:=0;
     end selfClusteringBondVolume;

     replaceable function selfClusteringBondHeatCapacityCp
        extends Modelica.Icons.Function;
            input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarHeatCapacity selfClusteringBondHeatCapacityCp;
     algorithm
       selfClusteringBondHeatCapacityCp:=0;
     end selfClusteringBondHeatCapacityCp;

      replaceable function specificAmountOfParticles
        "Amount of particles per mass of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
        output Real specificAmountOfSubstance(unit="mol/kg")
          "Amount of substance particles per its mass";
      algorithm
        specificAmountOfSubstance := 1/molarMassOfBaseMolecule(substanceData);
        annotation (Inline=true, smoothOrder=2);
      end specificAmountOfParticles;

      replaceable function specificAmountOfFreeBaseMolecule
        "Amount of substance free base molecule per mass of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
        output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
          "Amount of substance free base molecule per substance mass";
      algorithm
        specificAmountOfFreeBaseMolecule := 1/molarMassOfBaseMolecule(substanceData);
        annotation (Inline=true, smoothOrder=2);
      end specificAmountOfFreeBaseMolecule;


    /* replaceable function solution_temperature_
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
     extends Modelica.Icons.Function;
    input SubstanceData substanceData[:] "Data record of substances";
  input Modelica.Units.SI.MolarEnthalpy h
    "Molar enthalpy of solution (x*substances_h)";
  input Modelica.Units.SI.MoleFraction x[:]
    "Mole fractions of substances";
  input Modelica.Units.SI.Pressure p=100000 "Pressure";
  input Modelica.Units.SI.ElectricPotential v=0
    "Electric potential of the substance";
  input Modelica.Units.SI.MoleFraction I=0
    "Ionic strengh (mole fraction based)";

  output Modelica.Units.SI.Temperature T "Temperature";
    annotation (__Dymola_DymolaStoredErrors(thetext="/*replaceable function solution_temperature_
  \"Temperature of the solution from specific enthalpy and mass fractions of substances\"
    extends Modelica.Icons.Function;
   input SubstanceData substanceData[:] \"Data record of substances\";
 input Modelica.Units.SI.MolarEnthalpy h
   \"Molar enthalpy of solution (x*substances_h)\";
 input Modelica.Units.SI.MoleFraction x[:]
   \"Mole fractions of substances\";
 input Modelica.Units.SI.Pressure p=100000 \"Pressure\";
 input Modelica.Units.SI.ElectricPotential v=0
   \"Electric potential of the substance\";
 input Modelica.Units.SI.MoleFraction I=0
   \"Ionic strengh (mole fraction based)\";

 output Modelica.Units.SI.Temperature T \"Temperature\";
"));
end solution_temperature_;
*/

     replaceable function specificEnthalpy
       "Specific molar enthalpy of the substance with electric potential dependence"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
        "Specific enthalpy";

     algorithm

       specificEnthalpy := molarEnthalpy(
         substanceData,
         T,
         p,
         v,
         I)/
         molarMassOfBaseMolecule(substanceData);
     end specificEnthalpy;

     replaceable function specificVolume "Specific volume of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.SpecificVolume specificVolume "Specific volume";

     algorithm

      specificVolume := molarVolume(
           substanceData,
           T,
           p,
           v,
           I) /
         molarMassOfBaseMolecule(substanceData);
     end specificVolume;

      replaceable function specificHeatCapacityCp
      "Specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
        "Specific heat capacity at constant pressure";

      algorithm

      specificHeatCapacityCp := molarHeatCapacityCp(
           substanceData,
           T,
           p,
           v,
           I) /
         molarMassOfBaseMolecule(substanceData);
      end specificHeatCapacityCp;

     replaceable function temperature
      "Temperature of the substance from its enthalpy"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.Temperature T "Temperature";
     end temperature;

     replaceable function solution_temperature
      "Temperature of the solution from specific enthalpy and mass fractions of substances"
         extends Modelica.Icons.Function;
        input SubstanceData substanceData[:] "Data record of substances";
      input Modelica.Units.SI.SpecificEnthalpy h
        "Specific enthalpy of solution";
      input Modelica.Units.SI.MassFraction X[:]
        "Mass fractions of substances";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.Temperature T "Temperature";
     end solution_temperature;

     replaceable function density
          "Return density of the substance in the solution"
            extends Modelica.Icons.Function;
            input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.Density density "Density";
     end density;
      annotation (Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end StateOfMatter;

    package Incompressible "Incompressible as basic state of matter"
      extends StateOfMatter;

      redeclare record extends SubstanceData "Base substance data"

        parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=
             0.01801528 "Molar weight of the substance";

        parameter Modelica.Units.SI.ChargeNumberOfIon z=0
          "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

        parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")=
          DfG_25degC_1bar
          "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

        parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")=
          DfH_25degC
          "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

        parameter Modelica.Units.SI.ActivityCoefficient gamma=1
          "Activity coefficient of the substance";

        parameter Modelica.Units.SI.MolarHeatCapacity Cp=0
          "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
        parameter String References[1]={""}
          "References of these thermodynamical values";

        parameter Modelica.Units.SI.MolarEnergy DfG_25degC_1bar(displayUnit="kJ/mol")=
             0 "Obsolete parameter use DfH instead"
          annotation (Dialog(tab="Obsolete"));

        parameter Modelica.Units.SI.MolarEnergy DfH_25degC(displayUnit="kJ/mol")=
             0 "Obsolete parameter use DfG instead"
          annotation (Dialog(tab="Obsolete"));

        parameter Boolean SelfClustering=false
          "Pure substance is making clusters (weak bonds between molecules)";

        parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
          "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";
        //-20000
        parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
          "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

        parameter Modelica.Units.SI.Density density(displayUnit="kg/dm3") = 1000
          "Density of the pure substance (default density of water at 25degC)";

        //      parameter Modelica.SIunits.MolarHeatCapacity Cv = Cp
        //      "Molar heat capacity of the substance at constant volume";

        annotation (preferredView="info", Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end SubstanceData;

      redeclare function extends activityCoefficient
        "Return activity coefficient of the substance in the solution"
      algorithm
        activityCoefficient := substanceData.gamma;
      end activityCoefficient;

      redeclare function extends chargeNumberOfIon
        "Return charge number of the substance in the solution"
      algorithm
        chargeNumberOfIon := substanceData.z;
      end chargeNumberOfIon;

      redeclare function extends molarEnthalpyElectroneutral
        "Molar enthalpy of the pure electroneutral substance"
      algorithm
        //Molar enthalpy:
        // - temperature and pressure shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - dp*Vm)
        //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - (dp/dT)*Vm. As a result dH = dT*Cv - dp*Vm for incompressible substances.

        molarEnthalpyElectroneutral := substanceData.DfH + (T - 298.15)*
          substanceData.Cp;
        //   - (p - 100000) * molarVolumePure(substanceData,T,p,v,I);
      end molarEnthalpyElectroneutral;

      redeclare function extends molarEntropyPure
        "Molar entropy of the pure substance"
      algorithm
        //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
        //+ substanceData.Cv*log(T/298.15);

        //Molar entropy shift:
        // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
        // - pressure shift: with constant molar volume at constant temperature Vm*dP = -T*dS (small amount of work)
        molarEntropyPure := substanceData.Cp*log(T/298.15) - (molarVolumePure(
            substanceData,
            T,
            p,
            v,
            I)/T)*(p - 100000) + ((substanceData.DfH - substanceData.DfG)/298.15);

        //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
        //As data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
        //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
      end molarEntropyPure;

      redeclare function molarVolumePure
        "Molar volume of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
      algorithm
        molarVolumePure := substanceData.MolarWeight/substanceData.density;
        //incompressible
      end molarVolumePure;

      redeclare function extends molarHeatCapacityCp
        "Molar heat capacity of the substance at constant pressure"
      algorithm
        molarHeatCapacityCp := substanceData.Cp;
      end molarHeatCapacityCp;

      redeclare function extends molarMassOfBaseMolecule
        "Molar mass of the substance"
      algorithm
        molarMass := substanceData.MolarWeight;
      end molarMassOfBaseMolecule;

      redeclare function selfClustering
        "returns true if substance molecules are joining together to clusters"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        output Boolean selfClustering;
      algorithm
        selfClustering := substanceData.SelfClustering;
      end selfClustering;

      redeclare function selfClusteringBondEnthalpy
        "Enthalpy of joining two base molecules of the substance together to cluster"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
      algorithm
        selfClusteringEnthalpy := substanceData.SelfClustering_dH;
      end selfClusteringBondEnthalpy;

      redeclare function selfClusteringBondEntropy
        "Entropy of joining two base molecules of the substance together to cluster"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
      algorithm
        selfClusteringEntropy := substanceData.SelfClustering_dS;
      end selfClusteringBondEntropy;

      redeclare replaceable function specificAmountOfParticles
      "Amount of substance particles per its mass"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Real specificAmountOfSubstance(unit="mol/kg") "Amount of substance particles per its mass";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if not selfClustering(substanceData) then
          specificAmountOfSubstance := 1/substanceData.MolarWeight;
        else
          SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);

          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));

          specificAmountOfSubstance := 1/((SelfClustering_K + 1)*substanceData.MolarWeight);
        end if;
      end specificAmountOfParticles;

      redeclare function specificAmountOfFreeBaseMolecule
        "Amount of substance free base molecule per mass of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
        output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
          "Amount of substance free base molecule per substance mass";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if not selfClustering(substanceData) then
          specificAmountOfFreeBaseMolecule := 1/substanceData.MolarWeight;
        else
          SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);

          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
          specificAmountOfFreeBaseMolecule :=(1- SelfClustering_K/(SelfClustering_K + 1)) / substanceData.MolarWeight;
        end if;
        annotation (Inline=true, smoothOrder=2);
      end specificAmountOfFreeBaseMolecule;

      redeclare replaceable function specificEnthalpy
        "Specific molar enthalpy of the substance with electric potential dependence"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0
          "Ionic strengh (mole fraction based)";

        output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
          "Specific enthalpy";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if selfClustering(substanceData) then
          SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);
          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
        end if;

        specificEnthalpy := (molarEnthalpy(
            substanceData,
            T,
            p,
            v,
            I) + (if selfClustering(substanceData) then
          selfClusteringBondEnthalpy(substanceData)*SelfClustering_K/(
          SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);

        annotation (Inline=true, smoothOrder=2);
      end specificEnthalpy;

      redeclare replaceable function specificVolume
        "Specific volume of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0
          "Ionic strengh (mole fraction based)";

        output Modelica.Units.SI.SpecificVolume specificVolume
          "Specific volume";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if selfClustering(substanceData) then
          SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);
          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
        end if;

        specificVolume := (molarVolume(
            substanceData,
            T,
            p,
            v,
            I) + (if selfClustering(substanceData) then
          selfClusteringBondVolume(substanceData)*SelfClustering_K/(
          SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);
      end specificVolume;

      redeclare replaceable function specificHeatCapacityCp
        "Specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature T=298.15 "Temperature";
        input Modelica.Units.SI.Pressure p=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
          "Specific heat capacity at constant pressure";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if selfClustering(substanceData) then
          SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);
          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
        end if;

        specificHeatCapacityCp := (molarHeatCapacityCp(
            substanceData,
            T,
            p,
            v,
            I) + (if selfClustering(substanceData) then
          selfClusteringBondHeatCapacityCp(substanceData)*SelfClustering_K/(
          SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);

          //TODO: + selfClusteringBondEnthalpy * der(K/(K + 1))/der(T) .. if (selfClusteringBondHeatCapacityCp!=0)
      end specificHeatCapacityCp;

      redeclare function extends temperature
        "Temperature of substance from its enthalpy"
    protected
        Modelica.Units.SI.SpecificEnthalpy baseSpecificEnthalpy;
      algorithm

        baseSpecificEnthalpy := specificEnthalpy(
            substanceData,
            298.15,
            p,
            v,
            I);

        T := 298.15 + (h - baseSpecificEnthalpy)/specificHeatCapacityCp(
          substanceData);
      end temperature;

      redeclare function extends solution_temperature
        "Temperature of the solution from enthalpies os substances"
        // Modelica.Units.SI.MoleFraction x[size(X, 1)];
    protected
        Modelica.Units.SI.SpecificEnthalpy solution_h_base;
      /*  Modelica.Units.SI.SpecificHeatCapacity solution_Cp=sum(X[i]*
      substanceData[i].Cp/molarMassOfBaseMolecule(substanceData[i]) for
      i in 1:size(X, 1));*/
      algorithm
        solution_h_base := X*specificEnthalpy(
            substanceData,
            298.15,
            p,
            v,
            I);
        T := 298.15 + (h - solution_h_base)/(X*specificHeatCapacityCp(substanceData));
      end solution_temperature;

       redeclare function extends density
        "Return density of the substance in the solution"
       algorithm
        density := substanceData.density;
       end density;
      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Incompressible;
  end Interfaces;
  annotation (
preferredView="info",
version="1.4.0",
versionDate="2021-01-27",
dateModified = "2021-01-27 11:10:41Z",
conversion(
  from(version="1.3.1", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_1.4.mos",
        to="1.4.0"),
  from(version="1.3.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_1.4.mos",
        to="1.4.0"),
  from(version="1.2.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_1.4.mos",
        to="1.4.0"),
  from(version="1.1.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.1_to_1.4.mos",
        to="1.4.0"),
  from(version="1.0.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.0_to_1.4.mos",
        to="1.4.0")),
      uses( Modelica(version="4.0.0")),
  Documentation(revisions="<html>
<p>Copyright (c) 2021, Marek Matej&aacute;k, Ph.D. </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>", info="<html>
<p>During each electro-chemical process an <a href=\"modelica://Chemical.Components.Substance\">electro-chemical potential</a> of the substances is equilibrating and all thermodynamical properties of the homogenous chemical solutions are evaluated. </p>
<p>Processes: chemical reactions, gas dissolution, diffusion, membrane transports, osmotic fluxes, electrochemical cells, electrodes, ..</p>
<p>Please see the <a href=\"modelica://Chemical.UsersGuide.Overview\">overview</a>.</p>
</html>"));
end Chemical;

package ModelicaServices
  "ModelicaServices (Dymola implementation) - Models and functions used in the Modelica Standard Library requiring a tool specific implementation"
extends Modelica.Icons.Package;

package Machine

  final constant Real eps=1.e-15 "Biggest number such that 1.0 + eps = 1.0";

  final constant Real small=1.e-60
    "Smallest number such that small and -small are representable on the machine";

  final constant Real inf=1.e+60
    "Biggest Real number such that inf and -inf are representable on the machine";
  annotation (Documentation(info="<html>
<p>
Package in which processor specific constants are defined that are needed
by numerical algorithms. Typically these constants are not directly used,
but indirectly via the alias definition in
<a href=\"modelica://Modelica.Constants\">Modelica.Constants</a>.
</p>
</html>"));
end Machine;
annotation (
  Protection(access=Access.hide),
  preferredView="info",
  version="4.0.0",
  versionBuild=1,
  versionDate="2020-01-15",
  dateModified = "2023-01-25 10:08:00Z",
  uses(Modelica(version="4.0.0")),
  conversion(
    noneFromVersion="1.0",
    noneFromVersion="1.1",
    noneFromVersion="1.2",
    noneFromVersion="3.2.1",
 noneFromVersion="3.2.3"),
  Documentation(info="<html>
<p>
This package contains a set of functions and models to be used in the
Modelica Standard Library that requires a tool specific implementation.
These are:
</p>

<ul>
<li> <a href=\"modelica://ModelicaServices.Animation.Shape\">Animation.Shape</a>
     provides a 3-dim. visualization of elementary
     mechanical objects. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Animation.Surface\">Animation.Surface</a>
     provides a 3-dim. visualization of
     moveable parameterized surface. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Animation.Vector\">Animation.Vector</a>
     provides a 3-dim. visualization of a vector objects. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Vector\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Vector</a>
     via inheritance.</li>
	 
<li> <a href=\"modelica://ModelicaServices.ExternalReferences.loadResource\">ExternalReferences.loadResource</a>
     provides a function to return the absolute path name of an URI or a local file name. It is used in
<a href=\"modelica://Modelica.Utilities.Files.loadResource\">Modelica.Utilities.Files.loadResource</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Machine\">Machine</a>
     provides a package of machine constants. It is used in
<a href=\"modelica://Modelica.Constants\">Modelica.Constants</a>.</li>

<li> <a href=\"modelica://ModelicaServices.System.exit\">System.exit</a> provides a function to terminate the execution of the Modelica environment. It is used in <a href=\"modelica://Modelica.Utilities.System.exit\">Modelica.Utilities.System.exit</a> via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Types.SolverMethod\">Types.SolverMethod</a>
     provides a string defining the integration method to solve differential equations in
     a clocked discretized continuous-time partition (see Modelica 3.3 language specification).
     It is not yet used in the Modelica Standard Library, but in the Modelica_Synchronous library
     that provides convenience blocks for the clock operators of Modelica version &ge; 3.3.</li>
</ul>

<p>
This is the Dymola implementation.
</p>

<p>
Original version
<strong>Licensed by the Modelica Association under the 3-Clause BSD License</strong><br>
Copyright &copy; 2009-2020, Modelica Association and contributors.
</p>
<p>
<strong>Modifications licensed by Dassault Syst&egrave;mes, copyright &copy; 2009-2020.</strong>
</p>

</html>"));
end ModelicaServices;

package Physiolibrary
  "System biology, integrative physiology and pathophysiology modelling library"
  extends Modelica.Icons.Package;

  package Media
    extends Modelica.Icons.Package;

    package Water "Incompressible water with constant heat capacity"
      extends Interfaces.PartialMedium(
        ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
        final mediumName="Water",
        substanceNames={"H2O"},
        final singleState=true,
        final reducedX=true,
        final fixedX=false,
        reference_T=310.15,
        reference_p=101325,
        reference_X={1},
        SpecificEnthalpy(nominal=1.0e5),
        Density(start=1e3, nominal=1e3),
        AbsolutePressure(start=1.0e5, nominal=1.0e5),
        Temperature(
          min=273,
          max=350,
          start=310.15));


    protected
      package stateOfMatter = Chemical.Interfaces.Incompressible
        "Substances model to translate data into substance properties";

    public
      redeclare connector extends SubstancesPort
         Chemical.Interfaces.SubstancePort_a H2O "Free water molecule (in pure water is only cca 1 mol/kg free water molecules, other cca 54.5 mols are bounded together by hydrogen bonds)";
         Chemical.Interfaces.SubstancePort_a H "Free hydrogen ion H+";
         Chemical.Interfaces.SubstancePort_a O2 "Free oxygen molecule";
         Chemical.Interfaces.SubstancePort_a H2 "Free hydrogen molecule";
         Chemical.Interfaces.SubstancePort_a OH "Free hydroxide molecule OH-";
         Modelica.Electrical.Analog.Interfaces.Pin cathode "Electric cathode";
         Modelica.Electrical.Analog.Interfaces.Pin anode "Electric anode";
      end SubstancesPort;

    redeclare replaceable model extends SubstancesDecomposition "Just because Modelica in today version cannot work properly with nested connectors"
      Chemical.Interfaces.SubstancePort_a H2O annotation (Placement(transformation(extent={{90,-110},{110,-90}})));
      Chemical.Interfaces.SubstancePort_a H annotation (Placement(transformation(extent={{90,-70},{110,-50}}), iconTransformation(extent={{90,-70},{110,-50}})));
      Chemical.Interfaces.SubstancePort_a O2 annotation (Placement(transformation(extent={{90,50},{110,70}})));
      Chemical.Interfaces.SubstancePort_a H2 annotation (Placement(transformation(extent={{90,90},{110,110}})));
      Chemical.Interfaces.SubstancePort_a OH annotation (Placement(transformation(extent={{92,-10},{112,10}})));
      Modelica.Electrical.Analog.Interfaces.PositivePin cathode annotation (Placement(transformation(extent={{-10,90},{10,110}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin anode annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
    equation
      connect(H2O, substances.H2O) annotation (Line(points={{100,-100},{56,-100},{56,-70},{-80,-70},{-80,0},{-100,0}},
                                                                                                          color={158,66,200}));
      connect(H, substances.H) annotation (Line(points={{100,-60},{88,-60},{88,-46},{-76,-46},{-76,0},{-100,0}},      color={158,66,200}));
      connect(O2, substances.O2) annotation (Line(points={{100,60},{-72,60},{-72,0},{-100,0}},      color={158,66,200}));
      connect(H2, substances.H2) annotation (Line(points={{100,100},{22,100},{22,80},{-76,80},{-76,0},{-100,0},{-100,0}},     color={158,66,200}));
      connect(OH, substances.OH) annotation (Line(points={{102,0},{-100,0}},                        color={158,66,200}));
      connect(cathode, substances.cathode) annotation (Line(points={{0,100},{-100,100},{-100,0},{-100,0}}, color={0,0,255}));
      connect(anode, substances.anode) annotation (Line(points={{0,-100},{-100,-100},{-100,0}},                     color={0,0,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SubstancesDecomposition;

    public
      redeclare replaceable model extends ChemicalSolution
      protected
          Real I = 0 "mole-fraction-based ionic strength";
          Real logH,logOH,logO2,logH2,eq;
      equation
        v=substances.cathode.v-substances.anode.v;
        0=substances.cathode.i+substances.anode.i;
        _i = substances.cathode.i;
        _i + (-1)*Modelica.Constants.F*eq = 0 "electric current is flow of electrons";

        T = stateOfMatter.solution_temperature(
            {Substances.Water},
            h,
            {1},
            p);

        substances.H.u + substances.OH.u = substances.H2O.u "H+ + OH- <-> H2O";
        2*substances.H.q + 2*eq + 0.5*substances.O2.q = substances.H2O.q "2H+ + 2e- + (1/2)O2 <-> H2O";
        substances.H2.q = 2*substances.H.q + 2*eq "H2 <-> 2H+ + 2e-";

        logH=logOH;


        substances.H2O.u = stateOfMatter.electroChemicalPotentialPure( Substances.Water, T, p, v, I);
        substances.H.u = stateOfMatter.electroChemicalPotentialPure( Substances.H, T, p, v, I) +
                         Modelica.Constants.R*T*logH;
        substances.O2.u = stateOfMatter.electroChemicalPotentialPure( Substances.O2, T, p, v, I) +
                         Modelica.Constants.R*T*logO2;
        substances.H2.u = stateOfMatter.electroChemicalPotentialPure( Substances.H2, T, p, v, I) +
                         Modelica.Constants.R*T*logH2;
        substances.OH.u = stateOfMatter.electroChemicalPotentialPure( Substances.OH, T, p, v, I) +
                         Modelica.Constants.R*T*logOH;

        substances.H2O.h_outflow = stateOfMatter.molarEnthalpy( Substances.Water, T, p, v, I);
        substances.H.h_outflow = stateOfMatter.molarEnthalpy( Substances.H, T, p, v, I);
        substances.O2.h_outflow = stateOfMatter.molarEnthalpy( Substances.O2, T, p, v, I);
        substances.H2.h_outflow = stateOfMatter.molarEnthalpy( Substances.H2, T, p, v, I);
        substances.OH.h_outflow = stateOfMatter.molarEnthalpy( Substances.OH, T, p, v, I);

        enthalpyFromSubstances =
         substances.H2O.q * actualStream(substances.H2O.h_outflow) +
         substances.H.q * actualStream(substances.H.h_outflow) +
         substances.O2.q * actualStream(substances.O2.h_outflow) +
         substances.H2.q * actualStream(substances.H2.h_outflow) +
         substances.OH.q * actualStream(substances.OH.h_outflow)
          "enthalpy from substances";


        massFlows = {substances.H2O.q * Substances.Water.MolarWeight +
         substances.H.q * Substances.H.MolarWeight +
         substances.O2.q * Substances.O2.MolarWeight +
         substances.H2.q * Substances.H2.MolarWeight +
         substances.OH.q * Substances.OH.MolarWeight}
          "mass change of water";
      end ChemicalSolution;

      replaceable function extends specificEnthalpies_TpvI
      algorithm
         specificEnthalpy:=stateOfMatter.specificEnthalpy(
            {Substances.Water},
            T,p,v,I);
      end specificEnthalpies_TpvI;

    public
      redeclare model extends BaseProperties(final standardOrderComponents=true)
        "Base properties of medium"

      equation
        d = 1000;
        h = X*stateOfMatter.specificEnthalpy(
            {Substances.Water},
            T=T,
            p=p);
        u = h - p/d;
        MM = 1/(X*stateOfMatter.specificAmountOfParticles({Substances.Water}));
        R_s = 8.3144/MM;
        state.p = p;
        state.T = T;

      end BaseProperties;

      redeclare replaceable record extends ThermodynamicState
        "A selection of variables that uniquely defines the thermodynamic state"
        extends Modelica.Icons.Record;
        AbsolutePressure p "Absolute pressure of medium";
        Temperature T "Temperature of medium";
        annotation (Documentation(info="<html>

</html>"));
      end ThermodynamicState;

      redeclare replaceable function extends setState_pTX
        "Return thermodynamic state as function of p, T and composition X or Xi"
      algorithm
        state.p := p;
        state.T := T;
      end setState_pTX;

      redeclare replaceable function extends setState_phX
        "Return thermodynamic state as function of p, h and composition X or Xi"
      algorithm
        state.p := p;
        state.T := stateOfMatter.solution_temperature(
            {Substances.Water},
            h,
            {1},
            p);
      end setState_phX;

      redeclare function extends specificEnthalpy "Return specific enthalpy"
      algorithm
        h := stateOfMatter.specificEnthalpy(
            Substances.Water,
            T=state.T,
            p=state.p);
      end specificEnthalpy;

      redeclare function extends specificHeatCapacityCp
        "Return specific heat capacity at constant pressure"
      algorithm
        cp := stateOfMatter.specificHeatCapacityCp(
            Substances.Water,
            T=state.T,
            p=state.p);
        annotation (Documentation(info="<html>

</html>"));
      end specificHeatCapacityCp;

      redeclare function extends density
      algorithm
        d := 1000;
      end density;

      redeclare function extends temperature
      algorithm
        T := state.T;
      end temperature;

      redeclare function extends pressure
      algorithm
        p := state.p;
      end pressure;


      annotation (Documentation(info="<html>
<p>
This package is a <strong>template</strong> for <strong>new medium</strong> models. For a new
medium model just make a copy of this package, remove the
\"partial\" keyword from the package and provide
the information that is requested in the comments of the
Modelica source.
</p>
</html>",   revisions="<html>
<p><i>2021</i></p>
<p>Marek Matejak, http://www.physiolib.com </p>
<p>All rights reserved. </p>
</html>"));
    end Water;

    package Interfaces

      partial package PartialMedium

      extends Modelica.Media.Interfaces.PartialMedium;

        replaceable connector SubstancesPort

        annotation (
            Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}, initialScale = 0.2), graphics={  Rectangle(
                  extent={{-20,2},{20,-2}},
                  lineColor={158,66,200},
                  lineThickness=0.5),                                                                                                                                                                                                      Polygon(points={{-80,50},
                      {80,50},{100,30},{80,-40},{60,-50},{-60,-50},{-80,-40},{-100,30},{-80,50}},                                                                                                                                                                                                        lineColor = {0, 0, 0}, fillColor = {158,66,200}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-65, 25}, {-55, 15}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-5, 25}, {5, 15}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{55, 25}, {65, 15}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-35, -15}, {-25, -25}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{25, -15}, {35, -25}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}, initialScale = 0.2), graphics={  Polygon(points = {{-40, 25}, {40, 25}, {50, 15}, {40, -20}, {30, -25}, {-30, -25}, {-40, -20}, {-50, 15}, {-40, 25}}, lineColor = {0, 0, 0}, fillColor = {158,66,200}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-32.5, 7.5}, {-27.5, 12.5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-2.5, 12.5}, {2.5, 7.5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{27.5, 12.5}, {32.5, 7.5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-17.5, -7.5}, {-12.5, -12.5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{12.5, -7.5}, {17.5, -12.5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 70}, {150, 40}}, lineColor = {0, 0, 0}, textString = "%name")}),
            Documentation(info = "<html>
        <p>
        This connector defines the \"substances port\" that
        is used for cross-membrane transports of selected free base chemical substances.
        </p>
        </html>"));
        end SubstancesPort;

        replaceable partial model SubstancesDecomposition "Just because Modelica in today version cannot work properly with nested connectors"
          SubstancesPort substances annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
        equation

        end SubstancesDecomposition;

        replaceable partial model ChemicalSolution
          "Adaptor between selected free base chemical substances and medium substances"
          outer Modelica.Fluid.System system "System wide properties";

          SubstancesPort substances "free base chemical substances";
          Physiolibrary.Types.RealIO.PressureInput p "pressure";
          Physiolibrary.Types.RealIO.SpecificEnthalpyInput h "specific enthalpy";
          Physiolibrary.Types.RealIO.MassFractionInput X[nS] "mass fractions of medium substances";
          Physiolibrary.Types.RealIO.ElectricCurrentInput _i "electric current from substances";

          Physiolibrary.Types.RealIO.MassFlowRateInput substanceMassFlowsFromStream[nS] "flow of medium substances";
          Physiolibrary.Types.RealIO.MassOutput substanceMasses[nS](nominal=SubstanceFlowNominal) "mass od medium substances";

          parameter Types.Mass startSubstanceMasses[nS]=fill(Modelica.Constants.small,nS) "Initial value of medium substance masses";

          Physiolibrary.Types.RealIO.MassFlowRateOutput massFlows[nS](nominal=SubstanceFlowNominal) "mass flows trough substancesPort";
          Physiolibrary.Types.RealIO.TemperatureOutput T "temperature";
          Physiolibrary.Types.RealIO.HeatFlowRateOutput enthalpyFromSubstances "enthalpy from substances";

          Physiolibrary.Types.RealIO.ElectricPotentialOutput v "electric potential";

          Real logm[nS] "natutal logarithm of medium substance masses (as state variables)";
        initial equation
          //substanceMasses = startSubstanceMasses;
          logm = log(startSubstanceMasses);
        equation
        /* 
  enthalpyFromSubstances = 
   substances.*.q * actualStream(substances.*.h_outflow) "enthalpy from substances";
 

  massFlows = substances.\(*\).q .* MMb[\1];
*/

          //The main accumulation equation is "der(substanceMasses)= substanceMassFlowsFromStream + massFlows"
          // However, the numerical solvers can handle it in form of log(m) much better. :-)
          der(logm) = ((substanceMassFlowsFromStream + massFlows)./substanceMasses) "accumulation of substances=exp(logm) [kg]";
          substanceMasses = exp(logm);
        end ChemicalSolution;

        function i "Find index of substance"
          input String searchName "Name of substance to find in substanceNames";
          output Integer index "Index of searchName in substanceNames";
        algorithm
            index := -1;
            for i in 1:nS loop
              if ( Modelica.Utilities.Strings.isEqual(substanceNames[i], searchName)) then
               index := i;
              end if;
            end for;
            assert(index > 0, "Substance '" + searchName + "' is not present between Substances in Medium\n"
               + "Check parameters and medium model.");
        end i;


        constant Modelica.Units.SI.MassFlowRate SubstanceFlowNominal[nS]=ones(nS) "Nominal of substance flow";
        constant Modelica.Units.SI.SpecificEnthalpy SpecificEnthalpyNominal=-1E6 "Nominal of specific enthalpy";

         replaceable function specificEnthalpies_TpvI
              "Specific enthalpies of medium substances"
               input Modelica.Units.SI.Temperature T=298.15 "Temperature";
               input Modelica.Units.SI.Pressure p=100000 "Pressure";
               input Modelica.Units.SI.ElectricPotential v=0
                "Electric potential of the substance";
               input Modelica.Units.SI.MoleFraction I=0
                "Ionic strengh (mole fraction based)";
               output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy[nS]
                 "Specific enthalpies of medium substances";

         end specificEnthalpies_TpvI;

         function temperatureError  "To find u as temperature, where temperatureError(u,p,X,h)->0"
         extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
           input Real p;
           input Real X[nS];
           input Real h;
      protected
            Real hs[nS];
         algorithm
           hs:=specificEnthalpies_TpvI(u,p);
           y:=h-sum(hs[i]*X[i] for i in 1:nS);
         end temperatureError;

         partial function GetConcentration
           input ThermodynamicState state;
           output Types.Concentration C;
         end GetConcentration;

         partial function GetMassConcentration
           input ThermodynamicState state;
           output Types.MassConcentration R;
         end GetMassConcentration;

         partial function GetFraction
           input ThermodynamicState state;
           output Types.Fraction F;
         end GetFraction;

         partial function GetActivity
           input ThermodynamicState state;
           output Real A;
         end GetActivity;

        annotation (Documentation(revisions="<html>
<p><i>2021</i></p>
<p>Marek Matejak, http://www.physiolib.com </p>
<p>All rights reserved. </p>
</html>"));
      end PartialMedium;
    end Interfaces;

    package Substances

       constant Chemical.Interfaces.Incompressible.SubstanceData O2=
            Chemical.Substances.Oxygen_aqueous();

        constant Chemical.Interfaces.Incompressible.SubstanceData H2=
            Chemical.Substances.Hydrogen_aqueous();

        constant Chemical.Interfaces.Incompressible.SubstanceData H=
            Chemical.Substances.Proton_aqueous();

        constant Chemical.Interfaces.Incompressible.SubstanceData OH=
            Chemical.Substances.Hydroxide_aqueous();

        constant Chemical.Interfaces.Incompressible.SubstanceData Water=
             Chemical.Substances.Water_liquid();
    end Substances;
  end Media;

  package Fluid "Physiological fluids with static and dynamic properties"
    extends Modelica.Icons.Package;

    package Components
      extends Modelica.Icons.Package;

      model Conductor "Hydraulic resistor, where conductance=1/resistance"
        extends Physiolibrary.Fluid.Interfaces.OnePort;
        extends Physiolibrary.Icons.HydraulicResistor;
        parameter Boolean useConductanceInput = false "=true, if external conductance value is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        parameter Physiolibrary.Types.HydraulicConductance Conductance = 0 "Hydraulic conductance if useConductanceInput=false" annotation (
          Dialog(enable = not useConductanceInput));
        Physiolibrary.Types.RealIO.HydraulicConductanceInput cond(start = Conductance) = c if useConductanceInput annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 60})));
      protected
        Physiolibrary.Types.HydraulicConductance c;
        constant Boolean GenerateConductanceConnection = true;
      equation
        if not useConductanceInput and GenerateConductanceConnection then
          c = Conductance;
        end if;
        volumeFlowRate = c * (q_in.p - q_out.p);
        annotation (
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-220, -40}, {200, -80}}, lineColor = {127, 0, 0}, fillColor = {58, 117, 175}, fillPattern = FillPattern.Solid, textString = "%name")}),
          Documentation(revisions = "<html>
<p><i>2017-2018</i></p>
<p>Marek Matejak, http://www.physiolib.com </p>
<p>All rights reserved. </p>
</html>", info = "<html>
<p>This hydraulic conductance (resistance) element contains two connector sides. No hydraulic medium volume is changing in this element during simulation. That means that sum of flow in both connector sides is zero. The flow through element is determined by <b>Ohm&apos;s law</b>. It is used conductance (=1/resistance) because it could be numerical zero better then infinity in resistance. </p>
</html>"));
      end Conductor;

      model Resistor
        extends Physiolibrary.Fluid.Components.Conductor(final Conductance = 1 / Resistance, final useConductanceInput = false, final GenerateConductanceConnection = false);
        parameter Physiolibrary.Types.HydraulicResistance Resistance = Modelica.Constants.inf "Hydraulic conductance if useConductanceInput=false" annotation (
          Dialog(enable = not useResistanceInput));
        parameter Boolean useResistanceInput = false "=true, if external resistance value is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        Physiolibrary.Types.RealIO.HydraulicResistanceInput resistance(start = Resistance) = 1 / c if useResistanceInput annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 60})));
      equation
        if not useResistanceInput then
          c = 1 / Resistance;
        end if;
      end Resistor;

      model ElasticVessel "Elastic compartment as chemical solution envelop"
        extends Physiolibrary.Icons.ElasticBalloon;
        extends Physiolibrary.Fluid.Interfaces.Accumulation(final pressure_start = p_initial);
        parameter String stateName=getInstanceName();
        parameter Types.HydraulicCompliance Compliance = 1e+3
        "Compliance e.g. TidalVolume/TidalPressureGradient if useComplianceInput=false"                                                       annotation (
          Dialog(enable = not useComplianceInput));
        parameter Types.Volume ZeroPressureVolume = 1e-11
        "Functional Residual Capacity. Maximal fluid volume, that does not generate pressure if useV0Input=false"                                                   annotation (
          Dialog(enable = not useV0Input));
        //default = 1e-5 ml
        parameter Types.Pressure ExternalPressure = if isExternalPressureAbsolute then system.p_ambient else 0
        "External pressure if useExternalPressureInput=false."                                                                                                        annotation (
          Dialog(enable = not useExternalPressureInput));
        parameter Types.Volume ResidualVolume = 1e-9
        "Residual volume. Or maximal fluid volume, which generate negative collapsing pressure in linear model"                                              annotation (
          Dialog(tab = "Advanced", group = "Pressure-Volume relationship"));
        Types.Volume excessVolume
        "Additional cavity volume (=fluid volume + internal space volume), that generate pressure";
        parameter Boolean useV0Input = false "=true, if zero-pressure-fluid_volume input is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        Types.RealIO.VolumeInput zeroPressureVolume(start = ZeroPressureVolume) = zpv
        if useV0Input                                                                                                                           annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {-80, 80}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin={-70,90})));
        parameter Boolean useComplianceInput = false "=true, if compliance input is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        Types.RealIO.HydraulicComplianceInput compliance( start = Compliance) = c
        if useComplianceInput                                                                           annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 80}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin={0,90})));
        parameter Boolean useExternalPressureInput = false "=true, if external pressure input is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        parameter Boolean isExternalPressureAbsolute = false "external pressure as absolute pressure? Relative to ambient otherwise." annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));

        Types.RealIO.PressureInput externalPressure(start = ExternalPressure) = ep
        if useExternalPressureInput                                                                            annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {80, 80}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin={70,90})));
        Types.RealIO.VolumeOutput fluidVolume= volume annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {116, -60}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {100, -80})));
        parameter Boolean useSigmoidCompliance = false "sigmoid compliance e.g. lungs" annotation (
          Evaluate = true,
          choices(checkBox = true),
          Dialog(tab = "Advanced", group = "Pressure-Volume relationship"));
        parameter Types.Volume VitalCapacity = 0.00493
        "Relative volume capacity if useSigmoidCompliance"                                                annotation (
          Dialog(enable = useSigmoidCompliance, tab = "Advanced", group = "Pressure-Volume relationship"));
        parameter Types.Volume BaseTidalVolume = 0.000543
        "Base value of tidal volume"                                                   annotation (
          Dialog(enable = useSigmoidCompliance, tab = "Advanced", group = "Pressure-Volume relationship"));
        Types.RealIO.VolumeInput internalSpace(
          start=InternalSpace)=is if useInternalSpaceInput
        "additional internal volume (e.g. another inserted compartment inside)"
          annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=180,
              origin={100,8}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={90,60})));
        parameter Boolean useInternalSpaceInput = false "=true, if internal space input is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        parameter Types.Volume InternalSpace = 0
        "Internal space if there is no pressure gradient"                                          annotation (
          Dialog(tab = "Advanced", group = "Pressure-Volume relationship"));
        Types.Pressure relative_pressure;


      protected
        constant Boolean GenerateComplianceConnection = true;
        parameter Types.Pressure p_initial = system.p_ambient;
        parameter Types.Volume BaseMeanVolume = ZeroPressureVolume + BaseTidalVolume / 2
        "Point of equality with linear presentation such as (FunctionalResidualCapacity + TidalVolume/2)";

        Types.Pressure d_sigmoid = (BaseMeanVolume - ResidualVolume) * (VitalCapacity - (BaseMeanVolume - ResidualVolume)) / (c * VitalCapacity);
        Types.Pressure c_sigmoid = (BaseMeanVolume - ZeroPressureVolume) / c + d_sigmoid * log(VitalCapacity / (BaseMeanVolume - ResidualVolume) - 1);

        Types.Volume zpv;
        Types.Pressure ep;
        Types.HydraulicCompliance c;
        Types.Volume is;

      equation
      //elastic compartment
        if not useV0Input then
          zpv = ZeroPressureVolume;
        end if;
        if not useComplianceInput and GenerateComplianceConnection then
          c = Compliance;
        end if;
        if not useExternalPressureInput then
          ep = ExternalPressure;
        end if;
        if not useInternalSpaceInput then
          is = InternalSpace;
        end if;
        excessVolume = max(0, volume + is - zpv + InternalSpace) - InternalSpace;
        relative_pressure = pressure - (if isExternalPressureAbsolute then ep else ep + system.p_ambient);



        pressure = (if not useSigmoidCompliance
        then
          smooth(0,
          if noEvent(volume > ResidualVolume)
             then
                excessVolume / c
             else   (-(if isExternalPressureAbsolute then ep-system.p_ambient else ep) / log(Modelica.Constants.eps)) * log(max(Modelica.Constants.eps, volume / ResidualVolume)))
        else
          (-d_sigmoid * log(VitalCapacity / (volume - ResidualVolume) - 1)) + c_sigmoid)
        + (if isExternalPressureAbsolute then ep else ep + system.p_ambient);


        assert(volume > Modelica.Constants.eps, "Attempt to reach negative volume!");

        annotation (
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-280, -104}, {280, -142}}, lineColor = {127, 0, 0}, fillColor = {58, 117, 175}, fillPattern = FillPattern.Solid, textString = "%name")}),
          Documentation(revisions = "<html>
<p>2020 by Marek Matejak, http://www.physiolib.com </p>
</html>", info = "<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
      end ElasticVessel;

      model ElasticVesselElastance
        extends Physiolibrary.Fluid.Components.ElasticVessel(final Compliance = 1 / Elastance, final useComplianceInput = false, final GenerateComplianceConnection = false);
        parameter Physiolibrary.Types.HydraulicElastance Elastance = 1 "Elastance if useComplianceInput=false" annotation (
          Dialog(enable = not useComplianceInput));
        Types.RealIO.HydraulicElastanceInput elastance(start = Elastance) = 1 / c if useElastanceInput annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 80}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {20, 90})));
        parameter Boolean useElastanceInput = false "=true, if elastance input is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
      equation
        if not useElastanceInput then
          c = 1 / Elastance;
        end if;
      end ElasticVesselElastance;

      model Inertia "Inertia of the volumetric flow"
        extends Physiolibrary.Fluid.Interfaces.OnePort;
        //(q_in(m_flow(start=massFlow_start)));
        extends Physiolibrary.Icons.Inertance;
        parameter Physiolibrary.Types.MassFlowRate massFlow_start = 0 "Mass flow start value" annotation (
          Dialog(group = "Initialization"));
        parameter Physiolibrary.Types.HydraulicInertance I "Inertance";
      initial equation
        q_in.m_flow = massFlow_start;
      equation
        I * der(q_in.m_flow) = q_in.p - q_out.p;
        annotation (
          Documentation(info = "<html>
<p>Inertance I of the simple tube could be calculated as I=ro*l/A, where ro is fuid density, l is tube length and A is tube cross-section area.</p>
</html>", revisions = "<html>
<p><i>2017-2018</i></p>
<p>Marek Matejak, http://www.physiolib.com </p>
</html>"),Icon(graphics={  Text(extent = {{-212, -58}, {208, -98}}, lineColor = {127, 0, 0}, fillColor = {58, 117, 175}, fillPattern = FillPattern.Solid, textString = "%name")}));
      end Inertia;

      model IdealValve
        extends Icons.IdealValve;
        extends Physiolibrary.Fluid.Interfaces.OnePort;
        Boolean open(start = true) "Switching state";
        Real passableVariable(start = 0, final unit = "1") "Auxiliary variable for actual position on the ideal diode characteristic";
        /*  = 0: knee point
              < 0: below knee point, diode locking
              > 0: above knee point, diode conducting */
        parameter Physiolibrary.Types.HydraulicConductance _Gon(final min = 0, displayUnit = "l/(mmHg.min)") = 1.2501026264094e-02 "Forward state-on conductance (open valve conductance)" annotation (
          Dialog(enable = not useLimitationInputs));
        //= the same as resistance 1e-5 mmHg/(l/min)
        parameter Physiolibrary.Types.HydraulicConductance _Goff(final min = 0, displayUnit = "l/(mmHg.min)") = 1.2501026264094e-12 "Backward state-off conductance (closed valve conductance)" annotation (
          Dialog(enable = not useLimitationInputs));
        //= 1e-5 (l/min)/mmHg
        parameter Physiolibrary.Types.Pressure Pknee(final min = 0) = 0 "Forward threshold pressure";
        parameter Boolean useLimitationInputs = false "=true, if Gon and Goff are from inputs" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        Physiolibrary.Types.RealIO.HydraulicConductanceInput Gon(start = _Gon) = gon if useLimitationInputs "open valve conductance = infinity for ideal case" annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {-60, 100})));
        Physiolibrary.Types.RealIO.HydraulicConductanceInput Goff(start = _Goff) = goff if useLimitationInputs "closed valve conductance = zero for ideal case" annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {60, 100})));
      protected
        constant Boolean GenerateConductanceConnection = true;
        Physiolibrary.Types.HydraulicConductance gon;
        Physiolibrary.Types.HydraulicConductance goff;
        constant Physiolibrary.Types.Pressure unitPressure = 1;
        constant Physiolibrary.Types.VolumeFlowRate unitFlow = 1;
      equation
        if not useLimitationInputs and GenerateConductanceConnection then
          gon = _Gon;
          goff = _Goff;
        end if;
        open = passableVariable > Modelica.Constants.eps;
        dp = passableVariable * unitFlow * (if open then 1 / gon else 1) + Pknee;
        volumeFlowRate = passableVariable * unitPressure * (if open then 1 else goff) + goff * Pknee;
        annotation (
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-188, -100}, {196, -130}}, lineColor = {127, 0, 0}, fillPattern = FillPattern.Sphere, fillColor = {255, 85, 85}, textString = "%name")}),
          Documentation(info = "<html>
<p>Ideal Valve allows a volumetric flow in one direction in case of pressure gradient is greater. </p>
</html>", revisions = "<html>
</html>"));
      end IdealValve;

      model IdealValveResistance
        extends Physiolibrary.Fluid.Components.IdealValve(final _Gon = 1 / _Ron, final _Goff = 1 / _Roff, final useLimitationInputs = false, GenerateConductanceConnection = false);
        parameter Physiolibrary.Types.HydraulicResistance _Ron(displayUnit = "(mmHg.min)/l") = 79.993432449 "forward state resistance" annotation (
          Dialog(enable = not useResistanceInputs));
        parameter Physiolibrary.Types.HydraulicResistance _Roff = 799934324490.0 "Backward state-off resistance (closed valve resistance)" annotation (
          Dialog(enable = not useResistanceInputs));
        parameter Boolean useResistanceInputs = false "=true, if Ron and Roff are from inputs" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        Physiolibrary.Types.RealIO.HydraulicResistanceInput Ron(start = _Ron) = 1 / gon if useResistanceInputs "open valve resistancece = zero for ideal case" annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {-60, 100})));
        Physiolibrary.Types.RealIO.HydraulicResistanceInput Roff(start = _Roff) = 1 / goff if useResistanceInputs "closed valve resistance = infinity for ideal case" annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {60, 100})));
      equation
        if not useResistanceInputs then
          gon = 1 / _Ron;
          goff = 1 / _Roff;
        end if;
      end IdealValveResistance;
      annotation (
        Documentation(info = "<html>
<p>Main components for physiological fluid modeling.</p>
</html>"));
    end Components;

    package Interfaces
      extends Modelica.Icons.InterfacesPackage;

      connector FluidPort = Modelica.Fluid.Interfaces.FluidPort(redeclare
          replaceable package Medium =
            Physiolibrary.Media.Water);

      connector FluidPort_a "Hydraulical inflow connector"
        extends FluidPort;
        annotation (
          defaultComponentName = "port_a",
          Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 110}, {150, 50}}, textString = "%name")}),
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {127, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid)}),
          Documentation(info = "<html>
<p>
Connector with one flow signal of type Real.
</p>
</html>", revisions = "<html>
<p><i>2017-2018</i></p>
<p>Marek Matejak, marek@matfyz.cz </p>
</html>"));
      end FluidPort_a;

      connector FluidPort_b "Hydraulical outflow connector"
        extends FluidPort;
        annotation (
          defaultComponentName = "port_b",
          Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, 30}, {30, -30}}, lineColor = {127, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 110}, {150, 50}}, textString = "%name")}),
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {127, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {127, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
      end FluidPort_b;

      connector FluidPorts_a "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
        extends FluidPort;
        annotation (
          defaultComponentName = "ports_a",
          Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-50, -200}, {50, 200}}, initialScale = 0.2), graphics={  Text(extent = {{-75, 130}, {75, 100}}, textString = "%name"), Rectangle(extent = {{25, -100}, {-25, 100}}, lineColor = {127, 0, 0}), Ellipse(extent = {{-25, 90}, {25, 40}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-25, 25}, {25, -25}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-25, -40}, {25, -90}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid)}),
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-50, -200}, {50, 200}}, initialScale = 0.2), graphics={  Rectangle(extent = {{50, -200}, {-50, 200}}, lineColor = {127, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5), Ellipse(extent = {{-50, 180}, {50, 80}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, 50}, {50, -50}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, -80}, {50, -180}}, lineColor = {0, 0, 0}, fillColor = {127, 0, 0}, fillPattern = FillPattern.Solid)}));
      end FluidPorts_a;

      partial model OnePort "Hydraulical OnePort"
        replaceable package Medium = Media.Water constrainedby
        Media.Interfaces.PartialMedium                                                        "Medium model" annotation (
           choicesAllMatching = true);
        outer Modelica.Fluid.System system "System wide properties";
        FluidPort_a q_in(redeclare package Medium = Medium) "Inflow" annotation (
          Placement(transformation(extent = {{-114, -14}, {-86, 14}})));
        FluidPort_b q_out(redeclare package Medium = Medium) "Outflow" annotation (
          Placement(transformation(extent = {{86, -14}, {114, 14}})));
        Physiolibrary.Types.MassFlowRate massFlowRate "Mass flow";
        Physiolibrary.Types.VolumeFlowRate volumeFlowRate "Volume flow";
        Physiolibrary.Types.Pressure dp "Pressure gradient";
        Modelica.Units.SI.Density density(start = Medium.density_pTX(system.p_ambient, system.T_ambient, Medium.reference_X));
        parameter Boolean EnthalpyNotUsed = false annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(tab = "Advanced", group = "Performance"));
      equation
        q_in.m_flow + q_out.m_flow = 0;
        massFlowRate = q_in.m_flow;
        dp = q_in.p - q_out.p;
        q_in.Xi_outflow = inStream(q_out.Xi_outflow);
        q_in.C_outflow = inStream(q_out.C_outflow);
        q_out.Xi_outflow = inStream(q_in.Xi_outflow);
        q_out.C_outflow = inStream(q_in.C_outflow);
        volumeFlowRate * density = massFlowRate;
        if EnthalpyNotUsed then
          q_in.h_outflow = Medium.specificEnthalpy_pTX(system.p_ambient, system.T_ambient, Medium.reference_X);
          q_out.h_outflow = Medium.specificEnthalpy_pTX(system.p_ambient, system.T_ambient, Medium.reference_X);
          density = Medium.density_pTX(q_in.p, system.T_ambient, Medium.reference_X);
        else
          q_in.h_outflow = inStream(q_out.h_outflow);
          q_out.h_outflow = inStream(q_in.h_outflow);
          // medium density
          density = if q_in.m_flow >= 0 then Medium.density_phX(q_in.p, inStream(q_in.h_outflow), inStream(q_in.Xi_outflow)) else Medium.density_phX(q_out.p, inStream(q_out.h_outflow), inStream(q_out.Xi_outflow));
        end if;
      end OnePort;

      partial model CompositionSetup "Initial substances composition setup"
        replaceable package Medium = Media.Water
            constrainedby Media.Interfaces.PartialMedium "Medium model" annotation (
           choicesAllMatching = true);

        outer Modelica.Fluid.System system "System wide properties";
        parameter Modelica.Units.SI.MassFraction massFractions_start[:] = Medium.reference_X "* Masses of all base molecules. If size is nS-1 then last value is 1-sum(others). If size is nS then all values are scaled to sum==1." annotation (
          Dialog(enable = not use_concentration_start, group = "Initialization of medium composition"));
        parameter Real extraConcentration_start[Medium.nC] = Medium.C_default "Extra substance amounts per kilogram of solution"
          annotation(Dialog(group = "Initialization of medium composition"));
        parameter Modelica.Units.SI.Temperature temperature_start = system.T_ambient "Initial temperature" annotation (
          Dialog(group = "Initialization"));
        parameter Modelica.Units.SI.Pressure pressure_start = system.p_ambient "Initial pressure" annotation (
          Dialog(group = "Initialization"));

      protected
        parameter Modelica.Units.SI.MassFraction x_mass_start[Medium.nS] =
          if Medium.nS < 2 then {1} else
                if size(massFractions_start, 1) == Medium.nS - 1 then
                 cat(1, massFractions_start, {1 - sum(massFractions_start)})
                elseif size(massFractions_start, 1) == Medium.nS then
                  massFractions_start
                else
                 ones(Medium.nS) "Initial mass fractions of substances";
        parameter Real C_start[Medium.nC] = extraConcentration_start "Extra substance amounts per kilogram of solution";
        annotation (
          Icon(coordinateSystem(preserveAspectRatio = false)),
          Diagram(coordinateSystem(preserveAspectRatio = false)));
      end CompositionSetup;

      partial model Accumulation
        extends Physiolibrary.Fluid.Interfaces.CompositionSetup;

        parameter Integer nPorts = 0 "Number of hydraulic ports" annotation (
          Evaluate = true,
          Dialog(connectorSizing = true, group = "Ports"));
        Interfaces.FluidPorts_a q_in[nPorts](redeclare package Medium = Medium, each h_outflow(nominal=Medium.SpecificEnthalpyNominal)) annotation (
          Placement(transformation(extent = {{-10, -28}, {10, 28}}), iconTransformation(extent = {{-7, -26}, {7, 26}}, rotation = 180, origin = {-1, 0})));
        parameter Boolean useSubstances = false "=true, if substance ports are used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        parameter Boolean onElectricGround = false "=true, if electric potencial is zero" annotation (
          Evaluate = true,
          choices(checkBox = true));
        //,Dialog(group="Conditional inputs"));
        Medium.SubstancesPort substances if useSubstances annotation (
          Placement(transformation(extent={{-120,-20},{-80,20}}),      iconTransformation(extent={{-120,
                  -20},{-80,20}})));

        Medium.ChemicalSolution chemicalSolution(
          startSubstanceMasses = m_start,
          p = pressure,
          h = enthalpy / mass,
          X = if not Medium.reducedX then massFractions else cat(1, massFractions, {1 - sum(massFractions)}),
          _i = i)  if useSubstances;                              //enthalpy / mass,

        parameter Boolean use_mass_start = false "Use mass_start, otherwise volume_start" annotation (
          Evaluate = true,
          choices(checkBox = true),
          Dialog(group = "Initialization"));
        parameter Physiolibrary.Types.Volume volume_start=0.001   "Total volume of solution start value" annotation (
          HideResult = use_mass_start,
          Dialog(enable = not use_mass_start, group = "Initialization"));
        parameter Physiolibrary.Types.Mass mass_start(displayUnit="kg")=1     "Total mass of solution start value" annotation (
          HideResult = not use_mass_start,
          Dialog(enable = use_mass_start, group = "Initialization"));


        parameter Boolean useThermalPort = false "Is thermal port pressent?" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));

        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(T = Medium.temperature_phX(pressure, enthalpy / mass, massFractions), Q_flow = heatFromEnvironment) if useThermalPort annotation (
          Placement(transformation(extent = {{-70, -90}, {-50, -70}}), iconTransformation(extent={{-70,
                -110},{-50,-90}})));

      protected
        parameter Physiolibrary.Types.Mass tm_start(displayUnit = "kg") = if use_mass_start then mass_start else volume_start * Medium.density_pTX(pressure_start, temperature_start, x_mass_start) "If both mass_start and volume_start are filled";

        parameter Modelica.Units.SI.Mass m_start[Medium.nS] = tm_start * x_mass_start[1:Medium.nS];
        Modelica.Units.SI.ElectricCurrent i;
      public
        Physiolibrary.Types.HeatFlowRate heatFromEnvironment;

        Physiolibrary.Types.Enthalpy enthalpy( start = m_start * Medium.specificEnthalpies_TpvI(temperature_start,pressure_start));

        Physiolibrary.Types.Mass mass(start = tm_start);
        Physiolibrary.Types.MassFraction massFractions[Medium.nXi];
        Physiolibrary.Types.MassFraction xx_mass[nPorts, Medium.nXi] "Substance mass fraction per fluid port";

        Real xC_mass[nPorts, Medium.nC] "Extra substance in 1 kg of solution per fluid port";
        Real extraSubstanceAmounts[Medium.nC](start = tm_start * C_start) "Current amount of extra substances";
        Real extraSubstanceConcentrations[Medium.nC](start = C_start) "Current anount per kg of extra substances";

        Physiolibrary.Types.Volume volume;
        Physiolibrary.Types.Density density;
      protected
        Physiolibrary.Types.Pressure pressure;
        Physiolibrary.Types.RealIO.HeatFlowRateOutput enthalpyFromSubstances "Enthalpy inflow in substances connectors [J/s]";
        Physiolibrary.Types.RealIO.MassFlowRateOutput massFlows[Medium.nS](nominal=Medium.SubstanceFlowNominal);
        Physiolibrary.Types.RealIO.ElectricPotentialOutput v;

        Physiolibrary.Types.RealIO.MassFlowRateOutput substanceMassFlowsFromStream[Medium.nS](nominal=Medium.SubstanceFlowNominal);
        Physiolibrary.Types.RealIO.MassInput substanceMasses[Medium.nS](nominal=Medium.SubstanceFlowNominal);

      initial equation
      //  assert(abs(1 - sum(x_mass_start)) < 1e-5, "Sum of x_mass_start must be 1. (Composition initialization failed)");
      /* assert(
  not ((compositionType == Physiolibrary.Fluid.Interfaces.CompositionType.Concentration) and (size(concentration_start,1)==Medium.nS-2) and (Medium.nS<2) or 
  (Medium.zb[Medium.nS - 1]==0)), "Initial electroneutral concentration composition is not supported with this medium (try to use mass fractions)!");
*/
      /*  assert(
  not ((compositionType == Physiolibrary.Fluid.Interfaces.CompositionType.Concentration) and (size(concentration_start,1)>=Medium.nS-2)),
  "Initial concentration composition must have at least 
  -2 values!");
  */
        if not useSubstances then
          substanceMasses = m_start;
        end if;
        if Medium.reducedX then
          mass = tm_start;
        end if;

        enthalpy = m_start * Medium.specificEnthalpies_TpvI(temperature_start,pressure_start,v);

      equation

        if onElectricGround then
          v = 0;
        else
          i = 0;
        end if;
        if not useThermalPort then
          heatFromEnvironment = 0;
        end if;
        if useSubstances then
          connect(substances, chemicalSolution.substances);
          connect(chemicalSolution.massFlows, massFlows);
          connect(chemicalSolution.enthalpyFromSubstances, enthalpyFromSubstances);
          connect(chemicalSolution.substanceMasses, substanceMasses);
          connect(chemicalSolution.substanceMassFlowsFromStream, substanceMassFlowsFromStream);
          connect(v, chemicalSolution.v);
        else
          der(substanceMasses) = substanceMassFlowsFromStream;

          massFlows = zeros(Medium.nS);

          enthalpyFromSubstances = 0;

          if not onElectricGround then
          //both electric variables set to zero
            v = 0;
          else
            i = 0;
          end if;
        end if;

        substanceMassFlowsFromStream =  (if not Medium.reducedX then q_in.m_flow*xx_mass else cat(1, q_in.m_flow*xx_mass, {q_in.m_flow*(ones(nPorts) - xx_mass*ones(Medium.nXi))}));


        der(extraSubstanceAmounts) = q_in.m_flow * xC_mass;


        mass = sum(substanceMasses);

        massFractions = substanceMasses[1:Medium.nXi] ./ mass;

        der(enthalpy) = q_in.m_flow * actualStream(q_in.h_outflow) + enthalpyFromSubstances + heatFromEnvironment;

        volume = mass / density;
        density = Medium.density_phX(pressure, enthalpy / mass, massFractions);

        extraSubstanceConcentrations = extraSubstanceAmounts ./ mass;
        for i in 1:nPorts loop
          xx_mass[i, :] = actualStream(q_in[i].Xi_outflow);
          xC_mass[i, :] = actualStream(q_in[i].C_outflow);
          q_in[i].p = pressure;
          q_in[i].h_outflow = enthalpy / mass;
          q_in[i].Xi_outflow = massFractions;
          q_in[i].C_outflow  = extraSubstanceConcentrations;
        end for;

        annotation (
          Icon(coordinateSystem(preserveAspectRatio = false)),
          Diagram(coordinateSystem(preserveAspectRatio = false)));
      end Accumulation;

      partial model PartialAbsoluteSensor "Partial component to model a sensor that measures a potential variable"
        replaceable package Medium = Physiolibrary.Media.Water constrainedby
        Physiolibrary.Media.Interfaces.PartialMedium                                                                      "Medium in the sensor" annotation (
           choicesAllMatching = true);
        Modelica.Fluid.Interfaces.FluidPort_a
                    port(redeclare package Medium = Medium, m_flow(min = 0)) annotation (
          Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      equation
        port.m_flow = 0;
        port.h_outflow = Medium.h_default;
        port.Xi_outflow = Medium.X_default[1:Medium.nXi];
        port.C_outflow = Medium.C_default;
        annotation (
          Documentation(info = "<html>
<p>
Partial component to model an <strong>absolute sensor</strong>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <code>PartialFlowSensor</code> instead.
as signal.
</p>
</html>"));
      end PartialAbsoluteSensor;
    end Interfaces;

    package Sensors
      extends Modelica.Icons.SensorsPackage;

      model PressureMeasure "Hydraulic pressure at port"
        extends Physiolibrary.Icons.PressureMeasure;
        extends Fluid.Interfaces.PartialAbsoluteSensor;

        outer Modelica.Fluid.System system "System wide properties";
        parameter Boolean GetAbsolutePressure = false "if false then output pressure is relative to ambient pressure" annotation (
          Evaluate = true,
          choices(checkBox = true));
        Physiolibrary.Types.RealIO.PressureOutput pressure "Pressure" annotation (
          Placement(transformation(extent = {{40, -60}, {80, -20}})));
      equation
        pressure =port.p  - (if GetAbsolutePressure then 0 else system.p_ambient);

        annotation (
          Documentation(revisions = "<html>
        <p><i>2009-2018</i></p>
        <p>Marek Matejak, marek@matfyz.cz </p>
        </html>"));
      end PressureMeasure;

      model Power "Power as pressure multiplied by volumetric flow between ports"
        extends Physiolibrary.Fluid.Interfaces.OnePort;
        extends Modelica.Icons.RoundSensor;
        Physiolibrary.Types.RealIO.PowerOutput power "Actual power" annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, -60}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {0, 120})));
      equation
        q_out.p = q_in.p;
        power = q_in.p * volumeFlowRate;
        annotation (
          Documentation(revisions = "<html>
        <p><i>2009-2018</i></p>
        <p>Marek Matejak, marek@matfyz.cz </p>
        </html>"),
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-25, -11}, {34, -70}}, lineColor = {0, 0, 0}, textString = "V'")}));
      end Power;

      model Sphygmomanometer "Systolic, diastolic and mean pressure measurement with latency of measurement time"
        extends Physiolibrary.Icons.PressureMeasure;
        extends Fluid.Interfaces.PartialAbsoluteSensor;
        parameter Physiolibrary.Types.Time MeasurementTime = 2 "Measurement time period";
        outer Modelica.Fluid.System system "System wide properties";
        Physiolibrary.Types.RealIO.PressureOutput diastolic "Diastolic pressure" annotation (
          Placement(transformation(extent = {{40, -60}, {80, -20}})));
        Types.RealIO.PressureOutput systolic "Systolic pressure" annotation (
          Placement(transformation(extent = {{40, 0}, {80, 40}})));
        Types.RealIO.PressureOutput mean "Mean pressure" annotation (
          Placement(transformation(extent = {{40, -30}, {80, 10}})));
      protected
        Boolean b;
        discrete Physiolibrary.Types.Time t0;
        discrete Physiolibrary.Types.Pressure systolicMeassured, diastolicMeassured, meanMeasured;
        Physiolibrary.Types.Pressure pressure, diastolicRunning, systolicRunning;
        Real pressureInt;
      equation
        diastolic = diastolicMeassured;
        systolic = systolicMeassured;
        mean = meanMeasured;
        der(diastolicRunning) = if pressure < diastolicRunning then min(0, der(pressure)) else 0;
        der(systolicRunning) = if pressure > systolicRunning then max(0, der(pressure)) else 0;
        der(pressureInt) = pressure;
        b = time - pre(t0) >= MeasurementTime;
        when {initial(), b} then
          t0 = time;
          diastolicMeassured = pre(diastolicRunning);
          systolicMeassured = pre(systolicRunning);
          meanMeasured = pre(pressureInt) / MeasurementTime;
          reinit(diastolicRunning, pressure);
          reinit(systolicRunning, pressure);
          reinit(pressureInt, 0);
        end when;
        pressure =port.p  - system.p_ambient;

        annotation (
          Documentation(revisions = "<html>
        <p><i>2009-2018</i></p>
        <p>Marek Matejak, marek@matfyz.cz </p>
        </html>"));
      end Sphygmomanometer;
    end Sensors;
    annotation (
      Documentation(info = "<html>
<p>The main usage of this fluid domain is modeling of the cardio-vascular, respiratory and lymhpatic system in human physiology. And because there are no extreme thermodynamic conditions, the system can be really simple &mdash;it is only necessary to model conditions for ideal gases, for incompressible liquids, at normal liquid temperatures and with relative pressure 5-20kPa. This boring thermodynamic state leads to the very simple blocks of resistance,  pressure, volumetric flow, inertia and finally the block of blood accumulation in elastic comparments.</p>
</html>"));
  end Fluid;

  package Icons "Icons for physiological models"
    extends Modelica.Icons.Package;

    class CardioVascular
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/csv.png")}));
    end CardioVascular;

    partial class ElasticBalloon
      annotation (
        Icon(graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/elastic_vessel.svg")}));
    end ElasticBalloon;

    partial class HydraulicResistor
      annotation (
        Icon(graphics={  Bitmap(extent = {{-120, -42}, {120, 44}}, fileName = "modelica://Physiolibrary/Resources/Icons/resistor.svg")}));
    end HydraulicResistor;

    class PressureMeasure
      annotation (
        Icon(graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/pressureMeassure.png")}));
    end PressureMeasure;

    class Inertance
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/inertia.svg")}));
    end Inertance;

    class IdealValve
      annotation (
        Icon(graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/ideal_valve.svg")}));
    end IdealValve;
    annotation (
      Documentation(revisions = ""));
  end Icons;

  package Types "Physiological units with nominals"
    extends Modelica.Icons.Package;

    package Constants
      extends Modelica.Icons.SourcesPackage;

      block FractionConst "Constant signal of type Fraction"
        parameter Types.Fraction k "Constant Fraction output value";
        RealIO.FractionOutput y "Fraction constant" annotation (
          Placement(transformation(extent = {{40, -10}, {60, 10}}), iconTransformation(extent = {{40, -10}, {60, 10}})));
      equation
        y = k;
        annotation (
          defaultComponentName = "fraction",
          Diagram(coordinateSystem(extent = {{-40, -40}, {40, 40}})),
          Icon(coordinateSystem(extent = {{-40, -40}, {40, 40}}, preserveAspectRatio = false), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, radius = 10, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -44}, {100, -64}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "%name"), Text(extent = {{-40, 10}, {40, -10}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "Const")}));
      end FractionConst;
    end Constants;

    package RealIO
      extends Modelica.Icons.Package;

      connector ElectricCurrentInput = input ElectricCurrent "input ElectricCurrent as connector" annotation (
        defaultComponentName = "electriccurrent",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type ElectricCurrent.
            </p>
            </html>"));

      connector HeatFlowRateOutput = output HeatFlowRate "output HeatFlowRate as connector" annotation (
        defaultComponentName = "heatflowrate",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector MassInput = input Mass "input Mass as connector" annotation (
        defaultComponentName = "mass",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type Mass.
            </p>
            </html>"));

      connector MassOutput = output Mass "output Mass as connector" annotation (
        defaultComponentName = "mass",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector MassFractionInput = input MassFraction "input Mass Fraction as connector" annotation (
        defaultComponentName = "massFraction",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type MassFraction.
            </p>
            </html>"));

      connector MassFlowRateInput = input MassFlowRate "input MassFlowRate as connector" annotation (
        defaultComponentName = "massflowrate",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type MassFlowRate.
            </p>
            </html>"));

      connector MassFlowRateOutput = output MassFlowRate "output MassFlowRate as connector" annotation (
        defaultComponentName = "massflowrate",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector PressureInput = input Pressure "input Pressure as connector" annotation (
        defaultComponentName = "pressure",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type Pressure.
            </p>
            </html>"));

      connector PressureOutput = output Pressure "output Pressure as connector" annotation (
        defaultComponentName = "pressure",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector VolumeInput = input Volume "input Volume as connector" annotation (
        defaultComponentName = "volume",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type Volume.
            </p>
            </html>"));

      connector VolumeOutput = output Volume "output Volume as connector" annotation (
        defaultComponentName = "volume",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector TemperatureOutput = output Temperature "output Temperature as connector" annotation (
        defaultComponentName = "temperature",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector TimeOutput = output Time "output Time as connector" annotation (
        defaultComponentName = "time",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector ElectricPotentialOutput = output ElectricPotential "output ElectricPotential as connector" annotation (
        defaultComponentName = "electricpotential",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector FractionOutput = output Fraction "output Fraction as connector" annotation (
        defaultComponentName = "fraction",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector FrequencyInput = input Frequency "input Frequency as connector" annotation (
        defaultComponentName = "frequency",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type Frequency.
            </p>
            </html>"));

      connector HydraulicConductanceInput = input HydraulicConductance "input HydraulicConductance as connector" annotation (
        defaultComponentName = "hydraulicconductance",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type HydraulicConductance.
            </p>
            </html>"));

      connector HydraulicResistanceInput = input HydraulicResistance "input HydraulicResistance as connector" annotation (
        defaultComponentName = "hydraulicResistance",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type HydraulicResistance.
            </p>
            </html>"));

      connector HydraulicComplianceInput = input HydraulicCompliance "input HydraulicCompliance as connector" annotation (
        defaultComponentName = "hydrauliccompliance",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type HydraulicCompliance.
            </p>
            </html>"));

      connector HydraulicElastanceInput = input HydraulicElastance "input HydraulicElastance as connector" annotation (
        defaultComponentName = "hydraulicelastance",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type HydraulicElastance.
            </p>
            </html>"));

      connector HydraulicElastanceOutput = output HydraulicElastance "output HydraulicElastance as connector" annotation (
        defaultComponentName = "hydraulicelastance",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
          </p>
          </html>"));

      connector SpecificEnthalpyInput = input SpecificEnthalpy "input SpecificEnthalpy as connector" annotation (
        defaultComponentName = "specificEnthalpy",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type SpecificEnthalpy.
            </p>
            </html>"));

      connector PowerOutput = output Power "output Power as connector" annotation (
        defaultComponentName = "power",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Power.
          </p>
          </html>"));
    end RealIO;

    type Time = Modelica.Units.SI.Time(displayUnit = "min", nominal = 60);

    type Frequency = Modelica.Units.SI.Frequency(displayUnit = "1/min");

    type Mass = Modelica.Units.SI.Mass(displayUnit = "g", nominal = 1e-3, min = 0, max = Modelica.Constants.inf);

    type MassFraction = Modelica.Units.SI.MassFraction(nominal = 0.1, min = ModelicaServices.Machine.small, max = Modelica.Constants.inf);

    type MassFlowRate = Modelica.Units.SI.MassFlowRate(displayUnit = "mg/min", nominal = 0.001);

    type Density = Modelica.Units.SI.Density(displayUnit = "kg/l", nominal = 1e-3);

    type Pressure = Modelica.Units.SI.Pressure(displayUnit = "mmHg", nominal = 1e5);

    type Volume = Modelica.Units.SI.Volume(displayUnit = "ml", nominal = 1e-6, min = 0, max = Modelica.Constants.inf);

    type VolumeFlowRate = Modelica.Units.SI.VolumeFlowRate(displayUnit = "ml/min", nominal = 1e-6 / 60);

    replaceable type Concentration = Modelica.Units.SI.Concentration(displayUnit = "mmol/l", min = ModelicaServices.Machine.small, max = Modelica.Constants.inf) constrainedby Real;

    type MassConcentration = Modelica.Units.SI.MassConcentration(displayUnit = "mg/l", nominal = 1e-3, min = ModelicaServices.Machine.small, max = Modelica.Constants.inf);

    type Temperature = Modelica.Units.SI.Temperature(displayUnit = "degC", nominal = 1, min = 0);

    type HeatFlowRate = Modelica.Units.SI.HeatFlowRate(displayUnit = "kcal/min", nominal = 4186.8 / 60);

    type Power = Modelica.Units.SI.Power(displayUnit = "kcal/min", nominal = 4186.8 / 60);

    type SpecificEnthalpy = Modelica.Units.SI.SpecificEnthalpy(displayUnit = "kcal/kg", nominal = 1e5);

    type ElectricPotential = Modelica.Units.SI.ElectricPotential(displayUnit = "mV", nominal = 1e-3);

    type ElectricCurrent = Modelica.Units.SI.ElectricCurrent(displayUnit = "meq/min", nominal = 9.64853399 * 10 ^ 4 / 1000 / 60);

    type Fraction = Real(final quantity = "Fraction", final unit = "1", displayUnit = "%", nominal = 1e-2);

    type HydraulicConductance = Real(final quantity = "HydraulicConductance", final unit = "m3/(Pa.s)", displayUnit = "l/(mmHg.min)", nominal = 1e-3 / (133.322387415 * 60), min = 0);

    type HydraulicResistance = Real(final quantity = "HydraulicConductance", final unit = "(Pa.s)/m3", displayUnit = "(mmHg.min)/l", nominal = 1e+3 * 133.322387415 * 60, min = 0);

    type HydraulicCompliance = Real(final quantity = "HydraulicCompliance", final unit = "m3/Pa", displayUnit = "ml/mmHg", nominal = 1e-6 / 133.322387415);

    type HydraulicElastance = Real(final quantity = "HydraulicElastance", final unit = "Pa/m3", displayUnit = "mmHg/ml", nominal = 133.322387415 / 1e-6);

    type HydraulicInertance = Real(final quantity = "HydraulicInertance", final unit = "Pa.s2/kg", displayUnit = "mmHg.min2/g", nominal = 133.322387415 * 60 ^ 2 / 1e-3);

    type Enthalpy = Modelica.Units.SI.Enthalpy(displayUnit = "kcal", nominal = 4186.8) "Heat energy";
    annotation (
      Documentation(revisions = "<html>
        <p>Copyright (c) 2017-2018, Marek Matej&aacute;k, http://www.physiolib.com </p>
        </html>"));
  end Types;
  annotation (
    preferredView = "info",
    version = "3.0.0-beta1",
    versionDate = "2022-08-24",
    dateModified = "2022-08-24 17:14:41Z",
    uses(Modelica(version = "4.0.0"), Complex(version = "4.0.0"), Chemical(version = "1.4.0")),
    conversion(
      from(version = "BioChem-1.0.1", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertBioChem_1.0.1_to_Physiolibrary_2.3.mos", to = "3.0.0"),
      from(version = "0.4980", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_0.4980_to_2.3.mos", to = "3.0.0"),
      from(version = "1.0", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_1.0_to_3.0.mos", to = "3.0.0"),
      from(version = "1.1", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_1.1_to_3.0.mos", to = "3.0.0"),
      from(version = "1.2", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_1.2_to_3.0.mos", to = "3.0.0"),
      from(version = "2.0", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.0_to_3.0.mos", to = "3.0.0"),
      from(version = "2.1", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.1_to_3.0.mos", to = "3.0.0"),
      from(version = "2.1.0", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.1_to_3.0.mos", to = "3.0.0"),
      from(version = "2.1.1", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.1_to_3.0.mos", to = "3.0.0"),
      from(version = "2.1.2", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.1_to_3.0.mos", to = "3.0.0"),
      from(version = "2.2.0", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.2_to_3.0.mos", to = "3.0.0"),
      from(version = "2.3.0", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.3_to_3.0.mos", to = "3.0.0"),
      from(version = "2.3.1", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.3_to_3.0.mos", to = "3.0.0"),
      from(version = "2.3.2", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_2.3_to_3.0.mos", to = "3.0.0"),
      from(version = "3.0.0-alpha11", script = "modelica://Physiolibrary/Resources/Scripts/Dymola/ConvertPhysiolibrary_from_3.0.0-alpha11_to_3.0.0.mos", to = "3.0.0")),
    Documentation(revisions = "<html>
<p>Copyright (c) 2023, Marek Matej&aacute;k, Ph.D. </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>", info = "<html>
<ul>
<li>Web pages: <a href=\"http://www.physiolibrary.org\">www.physiolibrary.org</a></li>
<li><a href=\"modelica://Physiolibrary.UsersGuide.Overview\">Overview</a></li>
<li><a href=\"modelica://Physiolibrary.UsersGuide.Connectors\">Connectors</a></li>
<li><a href=\"modelica://Physiolibrary.UsersGuide.Contact\">Contact</a> </li>
</ul>
<p><br>The origin of this Modelica Physiolibrary was in the first version of our HumMod Golem Edition model implementation, where it was called HumMod.Library. As the successors of Guyton&apos;s Medical Physiology School write, the original HumMod model is &ldquo;The best, most complete, mathematical model of human physiology ever created&rdquo;.</p>
<p>We are also developing many types of smaller physiological models for use in medical education, so it was essential to separate this library from our HumMod Modelica implementation. This separation improves the quality of the next HumMod release and provides a useful Modelica library to modelers in this bioscience.</p>
<p>The library contains only carefully-chosen elementary physiological laws, which are the basis of more complex physiological processes.</p>
<p><br>Physiology is a very progressive discipline, that examines how the living body works. And it is no surprise that all processes in the human body are driven by physical laws of nature. The great challenge is to marry old empirical experiments with the &ldquo;new&rdquo; physical principles. Many teams and projects in the word deal with this formalization of physiology, for example: Physiome, SBML, EuroPhysiome, VPH, CellML etc. It is our hope that this library helps this unflagging effort of physiologists to exactly describe the processes.</p>
</html>"));
end Physiolibrary;

package Modelica "Modelica Standard Library - Version 4.0.0"
extends Modelica.Icons.Package;

  package Blocks
  "Library of basic input/output control blocks (continuous, discrete, logical, table blocks)"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;

    package Interfaces
    "Library of connectors and partial models for input/output blocks"
      extends Modelica.Icons.InterfacesPackage;

      connector RealInput = input Real "'input Real' as connector" annotation (
        defaultComponentName="u",
        Icon(graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            points={{-100.0,100.0},{100.0,0.0},{-100.0,-100.0}})},
          coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}},
            preserveAspectRatio=true,
            initialScale=0.2)),
        Diagram(
          coordinateSystem(preserveAspectRatio=true,
            initialScale=0.2,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            points={{0.0,50.0},{100.0,0.0},{0.0,-50.0},{0.0,50.0}}),
          Text(
            textColor={0,0,127},
            extent={{-10.0,60.0},{-10.0,85.0}},
            textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one input signal of type Real.
</p>
</html>"));

      connector RealOutput = output Real "'output Real' as connector" annotation (
        defaultComponentName="y",
        Icon(
          coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            points={{-100.0,100.0},{100.0,0.0},{-100.0,-100.0}})}),
        Diagram(
          coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            points={{-100.0,50.0},{0.0,0.0},{-100.0,-50.0}}),
          Text(
            textColor={0,0,127},
            extent={{30.0,60.0},{30.0,110.0}},
            textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one output signal of type Real.
</p>
</html>"));

      connector BooleanInput = input Boolean "'input Boolean' as connector"
        annotation (
        defaultComponentName="u",
        Icon(graphics={Polygon(
              points={{-100,100},{100,0},{-100,-100},{-100,100}},
              lineColor={255,0,255},
              fillColor={255,0,255},
              fillPattern=FillPattern.Solid)}, coordinateSystem(
            extent={{-100,-100},{100,100}},
            preserveAspectRatio=true,
            initialScale=0.2)),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            initialScale=0.2,
            extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{0,50},{100,0},{0,-50},{0,50}},
              lineColor={255,0,255},
              fillColor={255,0,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-10,85},{-10,60}},
              textColor={255,0,255},
              textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one input signal of type Boolean.
</p>
</html>"));

      connector BooleanOutput = output Boolean "'output Boolean' as connector"
        annotation (
        defaultComponentName="y",
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{-100,100},{100,0},{-100,-100},{-100,100}},
              lineColor={255,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{-100,50},{0,0},{-100,-50},{-100,50}},
              lineColor={255,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{30,110},{30,60}},
              textColor={255,0,255},
              textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one output signal of type Boolean.
</p>
</html>"));

      partial block SO "Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        RealOutput y "Connector of Real output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>
Block has one continuous Real output signal.
</p>
</html>"));

      end SO;

      partial block SISO "Single Input Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        RealInput u "Connector of Real input signal" annotation (Placement(
              transformation(extent={{-140,-20},{-100,20}})));
        RealOutput y "Connector of Real output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>
Block has one continuous Real input and one continuous Real output signal.
</p>
</html>"));
      end SISO;

      partial block SI2SO
        "2 Single Input / 1 Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        RealInput u1 "Connector of Real input signal 1" annotation (Placement(
              transformation(extent={{-140,40},{-100,80}})));
        RealInput u2 "Connector of Real input signal 2" annotation (Placement(
              transformation(extent={{-140,-80},{-100,-40}})));
        RealOutput y "Connector of Real output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));

        annotation (Documentation(info="<html>
<p>
Block has two continuous Real input signals u1 and u2 and one
continuous Real output signal y.
</p>
</html>"));

      end SI2SO;

      partial block SIMO "Single Input Multiple Output continuous control block"
        extends Modelica.Blocks.Icons.Block;
        parameter Integer nout=1 "Number of outputs";
        RealInput u "Connector of Real input signal" annotation (Placement(
              transformation(extent={{-140,-20},{-100,20}})));
        RealOutput y[nout] "Connector of Real output signals" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));

        annotation (Documentation(info="<html>
<p> Block has one continuous Real input signal and a
    vector of continuous Real output signals.</p>

</html>"));
      end SIMO;

      partial block partialBooleanSource
        "Partial source block (has 1 output Boolean signal and an appropriate default icon)"
        extends Modelica.Blocks.Icons.PartialBooleanBlock;

        Blocks.Interfaces.BooleanOutput y "Connector of Boolean output signal"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));

        annotation (
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                  100}}), graphics={
              Polygon(
                points={{-80,88},{-88,66},{-72,66},{-80,88}},
                lineColor={255,0,255},
                fillColor={255,0,255},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,66},{-80,-82}}, color={255,0,255}),
              Line(points={{-90,-70},{72,-70}}, color={255,0,255}),
              Polygon(
                points={{90,-70},{68,-62},{68,-78},{90,-70}},
                lineColor={255,0,255},
                fillColor={255,0,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{71,7},{85,-7}},
                lineColor=DynamicSelect({235,235,235}, if y then {0,255,0} else {235,235,235}),
                fillColor=DynamicSelect({235,235,235}, if y then {0,255,0} else {235,235,235}),
                fillPattern=FillPattern.Solid)}),
          Documentation(info="<html>
<p>
Basic block for Boolean sources of package Blocks.Sources.
This component has one continuous Boolean output signal y
and a 3D icon (e.g., used in Blocks.Logical library).
</p>
</html>"));

      end partialBooleanSource;

      partial block PartialNoise "Partial noise generator"
        import generator = Modelica.Math.Random.Generators.Xorshift128plus;
        import Modelica.Math.Random.Utilities.automaticLocalSeed;
        extends Modelica.Blocks.Interfaces.SO;

        // Main dialog menu
        parameter SI.Period samplePeriod(start=0.01)
          "Period for sampling the raw random numbers"
          annotation(Dialog(enable=enableNoise));

        // Advanced dialog menu: Noise generation
        parameter Boolean enableNoise = globalSeed.enableNoise
          "= true: y = noise, otherwise y = y_off"
          annotation(choices(checkBox=true),Dialog(tab="Advanced",group="Noise generation"));
        parameter Real y_off = 0.0
          "Sets y = y_off if enableNoise=false (or time<startTime, see below)"
          annotation(Dialog(tab="Advanced",group="Noise generation"));

        // Advanced dialog menu: Initialization
        parameter Boolean useGlobalSeed = true
          "= true: use global seed, otherwise ignore it"
          annotation(choices(checkBox=true),Dialog(tab="Advanced",group = "Initialization",enable=enableNoise));
        parameter Boolean useAutomaticLocalSeed = true
          "= true: use automatic local seed, otherwise use fixedLocalSeed"
          annotation(choices(checkBox=true),Dialog(tab="Advanced",group = "Initialization",enable=enableNoise));
        parameter Integer fixedLocalSeed = 1 "Local seed (any Integer number)"
          annotation(Dialog(tab="Advanced",group = "Initialization",enable=enableNoise and not useAutomaticLocalSeed));
        parameter SI.Time startTime = 0.0
          "Start time for sampling the raw random numbers"
          annotation(Dialog(tab="Advanced", group="Initialization",enable=enableNoise));
        final parameter Integer localSeed(fixed=false) "The actual localSeed";
      protected
        outer Modelica.Blocks.Noise.GlobalSeed globalSeed
          "Definition of global seed via inner/outer";
        parameter Integer actualGlobalSeed = if useGlobalSeed then globalSeed.seed else 0
          "The global seed, which is actually used";
        parameter Boolean generateNoise = enableNoise and globalSeed.enableNoise
          "= true, if noise shall be generated, otherwise no noise";

        // Declare state and random number variables
        Integer state[generator.nState] "Internal state of random number generator";
        discrete Real r "Random number according to the desired distribution";
        discrete Real r_raw "Uniform random number in the range (0,1]";

      initial equation
         localSeed = if useAutomaticLocalSeed then automaticLocalSeed(getInstanceName()) else fixedLocalSeed;
         pre(state) = generator.initialState(localSeed, actualGlobalSeed);
         r_raw = generator.random(pre(state));

      equation
        // Draw random number at sample times
        when generateNoise and sample(startTime, samplePeriod) then
          (r_raw, state) = generator.random(pre(state));
        end when;

        // Generate noise if requested
        y = if not generateNoise or time < startTime then y_off else r;

          annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics={
              Polygon(
                points={{-76,90},{-84,68},{-68,68},{-76,90}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-76,68},{-76,-80}}, color={192,192,192}),
              Line(points={{-86,-14},{72,-14}},
                                            color={192,192,192}),
              Polygon(
                points={{94,-14},{72,-6},{72,-22},{94,-14}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(visible = enableNoise,
                 points={{-76,-19},{-62,-19},{-62,-3},{-54,-3},{-54,-51},{-46,-51},{-46,
                    -29},{-38,-29},{-38,55},{-30,55},{-30,23},{-30,23},{-30,-37},{-20,
                    -37},{-20,-19},{-10,-19},{-10,-47},{0,-47},{0,35},{6,35},{6,49},{12,
                    49},{12,-7},{22,-7},{22,5},{28,5},{28,-25},{38,-25},{38,47},{48,47},
                    {48,13},{56,13},{56,-53},{66,-53}}),
              Text(
                extent={{-150,-110},{150,-150}},
                textString="%samplePeriod s"),
              Line(visible=not enableNoise,
                points={{-76,48},{72,48}}),
              Text(visible=not enableNoise,
                extent={{-75,42},{95,2}},
                textString="%y_off"),
              Text(visible=enableNoise and not useAutomaticLocalSeed,
                extent={{-92,20},{98,-22}},
                textColor={238,46,47},
                textString="%fixedLocalSeed")}),
          Documentation(info="<html>
<p>
Partial base class of noise generators defining the common features
of noise blocks.
</p>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end PartialNoise;
      annotation (Documentation(info="<html>
<p>
This package contains interface definitions for
<strong>continuous</strong> input/output blocks with Real,
Integer and Boolean signals. Furthermore, it contains
partial models for continuous and discrete blocks.
</p>

</html>",     revisions="<html>
<ul>
<li><em>June 28, 2019</em>
       by Thomas Beutlich:<br>
       Removed obsolete blocks.</li>
<li><em>Oct. 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       Added several new interfaces.</li>
<li><em>Oct. 24, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       RealInputSignal renamed to RealInput. RealOutputSignal renamed to
       output RealOutput. GraphBlock renamed to BlockIcon. SISOreal renamed to
       SISO. SOreal renamed to SO. I2SOreal renamed to M2SO.
       SignalGenerator renamed to SignalSource. Introduced the following
       new models: MIMO, MIMOs, SVcontrol, MVcontrol, DiscreteBlockIcon,
       DiscreteBlock, DiscreteSISO, DiscreteMIMO, DiscreteMIMOs,
       BooleanBlockIcon, BooleanSISO, BooleanSignalSource, MI2BooleanMOs.</li>
<li><em>June 30, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized a first version, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.</li>
</ul>
</html>"));
    end Interfaces;

    package Logical
    "Library of components with Boolean input and output signals"
      extends Modelica.Icons.Package;

      block Switch "Switch between two Real signals"
        extends Modelica.Blocks.Icons.PartialBooleanBlock;
        Blocks.Interfaces.RealInput u1 "Connector of first Real input signal"
          annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
        Blocks.Interfaces.BooleanInput u2 "Connector of Boolean input signal"
          annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
        Blocks.Interfaces.RealInput u3 "Connector of second Real input signal"
          annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
        Blocks.Interfaces.RealOutput y "Connector of Real output signal"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));

      equation
        y = if u2 then u1 else u3;
        annotation (
          defaultComponentName="switch1",
          Documentation(info="<html>
<p>The Logical.Switch switches, depending on the
logical connector u2 (the middle connector)
between the two possible input signals
u1 (upper connector) and u3 (lower connector).</p>
<p>If u2 is <strong>true</strong>, the output signal y is set equal to
u1, else it is set equal to u3.</p>
</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{12,0},{100,0}},
                color={0,0,127}),
              Line(points={{-100,0},{-40,0}},
                color={255,0,255}),
              Line(points={{-100,-80},{-40,-80},{-40,-80}},
                color={0,0,127}),
              Line(points={{-40,12},{-40,-12}},
                color={255,0,255}),
              Line(points={{-100,80},{-38,80}},
                color={0,0,127}),
              Line(points=DynamicSelect({{-38,80},{6,2}}, if u2 then {{-38,80},{6,2}} else {{-38,-80},{6,2}}),
                color={0,0,127},
                thickness=1),
              Ellipse(lineColor={0,0,255},
                pattern=LinePattern.None,
                fillPattern=FillPattern.Solid,
                extent={{2,-8},{18,8}})}));
      end Switch;
      annotation (Documentation(info="<html>
<p>
This package provides blocks with Boolean input and output signals
to describe logical networks. A typical example for a logical
network built with package Logical is shown in the next figure:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Logical/LogicalNetwork1.png\"
     alt=\"LogicalNetwork1.png\">
</p>

<p>
The actual value of Boolean input and/or output signals is displayed
in the respective block icon as \"circle\", where \"white\" color means
value <strong>false</strong> and \"green\" color means value <strong>true</strong>. These
values are visualized in a diagram animation.
</p>
</html>"),     Icon(graphics={Line(
              points={{-86,-22},{-50,-22},{-50,22},{48,22},{48,-22},{88,-24}},
              color={255,0,255})}));
    end Logical;

    package Math
    "Library of Real mathematical functions as input/output blocks"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.Package;

      block Feedback "Output difference between commanded and feedback input"

        Interfaces.RealInput u1 "Commanded input" annotation (Placement(transformation(extent={{-100,
                  -20},{-60,20}})));
        Interfaces.RealInput u2 "Feedback input" annotation (Placement(transformation(
              origin={0,-80},
              extent={{-20,-20},{20,20}},
              rotation=90)));
        Interfaces.RealOutput y annotation (Placement(transformation(extent={{80,-10},
                  {100,10}})));

      equation
        y = u1 - u2;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes output <strong>y</strong> as <em>difference</em> of the
commanded input <strong>u1</strong> and the feedback
input <strong>u2</strong>:
</p>
<blockquote><pre>
<strong>y</strong> = <strong>u1</strong> - <strong>u2</strong>;
</pre></blockquote>
<p>
Example:
</p>
<blockquote><pre>
   parameter:   n = 2

results in the following equations:

   y = u1 - u2
</pre></blockquote>

</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Ellipse(
                lineColor={0,0,127},
                fillColor={235,235,235},
                fillPattern=FillPattern.Solid,
                extent={{-20,-20},{20,20}}),
              Line(points={{-60,0},{-20,0}}, color={0,0,127}),
              Line(points={{20,0},{80,0}}, color={0,0,127}),
              Line(points={{0,-20},{0,-60}}, color={0,0,127}),
              Text(extent={{-14,-94},{82,0}}, textString="-"),
              Text(
                textColor={0,0,255},
                extent={{-150,40},{150,80}},
                textString="%name")}));
      end Feedback;

      block Add "Output the sum of the two inputs"
        extends Interfaces.SI2SO;

        parameter Real k1=+1 "Gain of input signal 1";
        parameter Real k2=+1 "Gain of input signal 2";

      equation
        y = k1*u1 + k2*u2;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes output <strong>y</strong> as <em>sum</em> of the
two input signals <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
<strong>y</strong> = k1*<strong>u1</strong> + k2*<strong>u2</strong>;
</pre></blockquote>
<p>
Example:
</p>
<blockquote><pre>
   parameter:   k1= +2, k2= -3

results in the following equations:

   y = 2 * u1 - 3 * u2
</pre></blockquote>

</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-100,60},{-74,24},{-44,24}}, color={0,0,127}),
              Line(points={{-100,-60},{-74,-24},{-44,-24}}, color={0,0,127}),
              Ellipse(lineColor={0,0,127}, extent={{-50,-50},{50,50}}),
              Line(points={{50,0},{100,0}}, color={0,0,127}),
              Text(extent={{-40,40},{40,-40}}, textString="+"),
              Text(extent={{-100,52},{5,92}}, textString="%k1"),
              Text(extent={{-100,-92},{5,-52}}, textString="%k2")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={         Line(points={{50,0},{100,0}},
                color={0,0,255}),                                        Line(
                points={{50,0},{100,0}}, color={0,0,127})}));
      end Add;

      block Mean "Calculate mean over period 1/f"
        extends Modelica.Blocks.Interfaces.SISO;
        parameter SI.Frequency f(start=50) "Base frequency";
        parameter Real x0=0 "Start value of integrator state";
        parameter Boolean yGreaterOrEqualZero=false
          "= true, if output y is guaranteed to be >= 0 for the exact solution"
          annotation (Evaluate=true, Dialog(tab="Advanced"));
      protected
        parameter SI.Time t0(fixed=false) "Start time of simulation";
        Real x "Integrator state";
        discrete Real y_last "Last sampled mean value";
      initial equation
        t0 = time;
        x = x0;
        y_last = 0;
      equation
        der(x) = u;
        when sample(t0 + 1/f, 1/f) then
          y_last = if not yGreaterOrEqualZero then f*pre(x) else max(0.0, f*pre(x));
          reinit(x, 0);
        end when;
        y = y_last;
        annotation (Documentation(info="<html>
<p>
This block calculates the mean of the input signal u over the given period 1/f:
</p>
<blockquote><pre>
1 T
- &int; u(t) dt
T 0
</pre></blockquote>
<p>
Note: The output is updated after each period defined by 1/f.
</p>

<p>
If parameter <strong>yGreaterOrEqualZero</strong> in the Advanced tab is <strong>true</strong> (default = <strong>false</strong>),
then the modeller provides the information that the mean of the input signal is guaranteed
to be &ge; 0 for the exact solution. However, due to inaccuracies in the numerical integration scheme,
the output might be slightly negative. If this parameter is set to true, then the output is
explicitly set to 0.0, if the mean value results in a negative value.
</p>
</html>"),     Icon(graphics={Text(
                extent={{-80,60},{80,20}},
                textString="mean"), Text(
                extent={{-80,-20},{80,-60}},
                textString="f=%f")}));
      end Mean;
      annotation (Documentation(info="<html>
<p>
This package contains basic <strong>mathematical operations</strong>,
such as summation and multiplication, and basic <strong>mathematical
functions</strong>, such as <strong>sqrt</strong> and <strong>sin</strong>, as
input/output blocks. All blocks of this library can be either
connected with continuous blocks or with sampled-data blocks.
</p>
</html>",     revisions="<html>
<ul>
<li><em>August 24, 2016</em>
       by Christian Kral: added WrapAngle</li>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       New blocks added: RealToInteger, IntegerToReal, Max, Min, Edge, BooleanChange, IntegerChange.</li>
<li><em>August 7, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized (partly based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist).
</li>
</ul>
</html>"),     Icon(graphics={Line(
              points={{-80,-2},{-68.7,32.2},{-61.5,51.1},{-55.1,64.4},{-49.4,72.6},
                  {-43.8,77.1},{-38.2,77.8},{-32.6,74.6},{-26.9,67.7},{-21.3,57.4},
                  {-14.9,42.1},{-6.83,19.2},{10.1,-32.8},{17.3,-52.2},{23.7,-66.2},
                  {29.3,-75.1},{35,-80.4},{40.6,-82},{46.2,-79.6},{51.9,-73.5},{
                  57.5,-63.9},{63.9,-49.2},{72,-26.8},{80,-2}},
              color={95,95,95},
              smooth=Smooth.Bezier)}));
    end Math;

    package Noise "Library of noise blocks"
      extends Modelica.Icons.Package;

      model GlobalSeed
        "Defines global settings for the blocks of sublibrary Noise, especially a global seed value is defined"
        parameter Boolean enableNoise = true
          "= true, if noise blocks generate noise as output; = false, if they generate a constant output"
          annotation(choices(checkBox=true));
        parameter Boolean useAutomaticSeed = false
          "= true, choose a seed by system time and process id; = false, use fixedSeed"
          annotation(choices(checkBox=true));
        parameter Integer fixedSeed = 67867967
          "Fixed global seed for random number generators (if useAutomaticSeed = false)"
            annotation(Dialog(enable=not useAutomaticSeed));
        final parameter Integer seed(fixed=false) "Actually used global seed";
        final parameter Integer id_impure(fixed=false)
          "ID for impure random number generators Modelica.Math.Random.Utilities.impureXXX" annotation(HideResult=true);
      initial equation
        seed = if useAutomaticSeed then Modelica.Math.Random.Utilities.automaticGlobalSeed() else fixedSeed;
        id_impure = Modelica.Math.Random.Utilities.initializeImpureRandom(seed);

        annotation (
          defaultComponentName="globalSeed",
          defaultComponentPrefixes="inner",
          missingInnerMessage="
Your model is using an outer \"globalSeed\" component but
an inner \"globalSeed\" component is not defined and therefore
a default inner \"globalSeed\" component is introduced by the tool.
To change the default setting, drag Noise.GlobalSeed
into your model and specify the seed.
",     Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                               graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
                                              Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              textColor={0,0,255}),
              Line(visible = enableNoise,
                   points={{-73,-15},{-59,-15},{-59,1},{-51,1},{-51,-47},{-43,-47},{-43,
                    -25},{-35,-25},{-35,59},{-27,59},{-27,27},{-27,27},{-27,-33},{-17,-33},{-17,-15},{-7,-15},{-7,-43},{3,
                    -43},{3,39},{9,39},{9,53},{15,53},{15,-3},{25,-3},{25,9},{31,9},{31,
                    -21},{41,-21},{41,51},{51,51},{51,17},{59,17},{59,-49},{69,-49}},
                  color={215,215,215}),
              Text(visible=enableNoise and not useAutomaticSeed,
                extent={{-90,-4},{88,-30}},
                textColor={255,0,0},
                textString="%fixedSeed"),
              Line(visible = not enableNoise,
                points={{-80,-4},{84,-4}},
                color={215,215,215}),
              Text(visible=enableNoise and not useAutomaticSeed,
                extent={{-84,34},{94,8}},
                textColor={255,0,0},
                textString="fixedSeed =")}),
          Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>",       info="<html>
<p>
When using one of the blocks of sublibrary <a href=\"modelica://Modelica.Blocks.Noise\">Noise</a>,
on the same or a higher hierarchical level, Noise.GlobalSeed
must be dragged resulting in a declaration
</p>

<blockquote><pre>
<strong>inner</strong> Modelica.Blocks.Noise.GlobalSeed globalSeed;
</pre></blockquote>

<p>
The GlobalSeed block provides global options for all Noise blocks of the same or a lower
hierarchical level. The following options can be selected:
</p>

<blockquote>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Icon</th>
    <th>Description</th></tr>

<tr><td> <img src=\"modelica://Modelica/Resources/Images/Blocks/Noise/GlobalSeed_FixedSeed.png\"> </td>
    <td> <strong>useAutomaticSeed=false</strong> (= default):<br>
         A fixed global seed is defined with Integer parameter fixedSeed. The value of fixedSeed
         is displayed in the icon. By default all Noise blocks use fixedSeed for initialization of their
         pseudo random number generators, in combination with a local seed defined for every instance
         separately. Therefore, whenever a simulation is performed with the
         same fixedSeed exactly the same noise is generated in all instances of the Noise
         blocks (provided the settings of these blocks are not changed as well).<br>
         This option can be used (a) to design a control system (e.g. by parameter optimization) and keep the same
         noise for all simulations, or (b) perform Monte Carlo Simulations where
         fixedSeed is changed from the environment for every simulation, in order to
         produce different noise at every simulation run.</td></tr>

<tr><td> <img src=\"modelica://Modelica/Resources/Images/Blocks/Noise/GlobalSeed_AutomaticSeed.png\"> </td>
    <td> <strong>useAutomaticSeed=true</strong>:<br>
         An automatic global seed is computed by using the ID of the process in which the
         simulation takes place and the current local time. As a result, the global seed
         is changed automatically for every new simulation, including parallelized
         simulation runs. This option can be used to perform Monte Carlo Simulations
         with minimal effort (just performing many simulation runs) where
         every simulation run uses a different noise.</td></tr>

<tr><td> <img src=\"modelica://Modelica/Resources/Images/Blocks/Noise/GlobalSeed_NoNoise.png\"> </td>
    <td> <strong>enableNoise=false</strong>:<br>
         The noise in all Noise instances is switched off and the blocks output a constant
         signal all the time (usually zero). This option is useful, if a model shall be
         tested without noise and the noise shall be quickly turned off or on.</td></tr>
</table>
</blockquote>

<p>
Additionally, the globalSeed instance calls function
<a href=\"modelica://Modelica.Math.Random.Utilities.initializeImpureRandom\">initializeImpureRandom</a>
to initialize the impure random number generators
(<a href=\"modelica://Modelica.Math.Random.Utilities.impureRandom\">impureRandom</a> and
<a href=\"modelica://Modelica.Math.Random.Utilities.impureRandomInteger\">impureRandomInteger</a>).
The return value of this function is stored in parameter <strong>id_impure</strong>. Whenever one of the impure
random number generators need to be called, \"globalSeed.id_impure\" has to be given as input argument.
</p>

<p>
Note, the usage of this block is demonstrated with examples
<a href=\"modelica://Modelica.Blocks.Examples.Noise.AutomaticSeed\">AutomaticSeed</a> and
<a href=\"modelica://Modelica.Blocks.Examples.Noise.ImpureGenerator\">ImpureGenerator</a>.
</p>

<p>
Please note that only one globalSeed instance may be defined in the model due to the initialization
of the impure random number generators with <a href=\"modelica://Modelica.Math.Random.Utilities.initializeImpureRandom\">initializeImpureRandom</a>!
So, the block will usually reside on the top level of the model.
</p>
</html>"));
      end GlobalSeed;

      block UniformNoise "Noise generator with uniform distribution"
        import distribution = Modelica.Math.Distributions.Uniform.quantile;
        extends Modelica.Blocks.Interfaces.PartialNoise;

        // Main dialog menu
        parameter Real y_min(start=0) "Lower limit of y" annotation(Dialog(enable=enableNoise));
        parameter Real y_max(start=1) "Upper limit of y" annotation(Dialog(enable=enableNoise));

      initial equation
         r = distribution(r_raw, y_min, y_max);

      equation
        // Draw random number at sample times
        when generateNoise and sample(startTime, samplePeriod) then
          r = distribution(r_raw, y_min, y_max);
        end when;

          annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics={
              Line(visible=enableNoise,
                points={{-76,60},{78,60}}, color={95,95,95},
                pattern=LinePattern.Dot),
              Line(visible=enableNoise,
                points={{-76,-60},{78,-60}}, color={95,95,95},
                pattern=LinePattern.Dot),
              Text(visible=enableNoise,
                extent={{-70,94},{95,64}},
                textColor={175,175,175},
                textString="%y_max"),
              Text(visible=enableNoise,
                extent={{-70,-64},{95,-94}},
                textColor={175,175,175},
                textString="%y_min")}),
          Documentation(info="<html>
<p>
A summary of the common properties of the noise blocks is provided in the documentation of package
<a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a>.
This UniformNoise block generates reproducible, random noise at its output according to a uniform distribution.
This means that random values are uniformly distributed within the range defined by parameters
y_min and y_max (see example <a href=\"modelica://Modelica.Blocks.Examples.Noise.UniformNoiseProperties\">Noise.UniformNoiseProperties</a>).
By default, two or more instances produce different, uncorrelated noise at the same time instant.
The block can only be used if on the same or a higher hierarchical level,
model <a href=\"modelica://Modelica.Blocks.Noise.GlobalSeed\">Blocks.Noise.GlobalSeed</a>
is dragged to provide global settings for all instances.
</p>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end UniformNoise;
      annotation (Icon(graphics={Line(
          points={{-84,0},{-54,0},{-54,40},{-24,40},{-24,-70},{6,-70},{6,80},
              {36,80},{36,-20},{66,-20},{66,60}})}), Documentation(info="<html>
<p>
This sublibrary contains blocks that generate <strong>reproducible noise</strong> with pseudo random
numbers. Reproducibility is important when designing control systems,
either manually or with optimization methods (for example when changing a parameter or a component
of a control system and re-simulating, it is important that the noise does not change, because
otherwise it is hard to determine whether the changed control system or the differently
computed noise has changed the behaviour of the controlled system).
Many examples how to use the Noise blocks are provided in sublibrary
<a href=\"modelica://Modelica.Blocks.Examples.Noise\">Blocks.Examples.Noise</a>.
</p>

<h4>Global Options</h4>

<p>
When using one of the blocks of this sublibrary, on the same or a higher level,
block <a href=\"modelica://Modelica.Blocks.Noise.GlobalSeed\">Noise.GlobalSeed</a>
must be dragged resulting in a declaration
</p>

<blockquote><pre>
<strong>inner</strong> Modelica.Blocks.Noise.GlobalSeed globalSeed;
</pre></blockquote>

<p>
This block is used to define global options that hold for all Noise block
instances (such as a global seed for initializing the random number generators,
and a flag to switch off noise). Furthermore, the impure random number generator
<a href=\"modelica://Modelica.Math.Random.Utilities.impureRandom\">impureRandom</a> is initialized here.
</p>

<p>
Please note that only one globalSeed instance may be defined in the model due to the initialization
of the impureRandom(..) random number generator! So, the block will usually reside on the top level of the model.
</p>

<h4>Parameters that need to be defined</h4>

<p>
When using a noise block of this package, at a minimum the following parameters must be defined:
</p>

<blockquote>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Parameter</th>
    <th>Description</th></tr>

<tr><td> samplePeriod </td>
    <td> Random values are drawn periodically at the sample rate in [s]
         defined with this parameter (time events are generated at the sample instants).
         Between sample instants, the output y is kept constant.</td></tr>

<tr><td> distribution data </td>
    <td> Every noise block in this package needs additional data to describe the respective
         distribution. A random number distribution maps the drawn random numbers
         from the range 0.0 ... 1.0, to the desired range and distribution.
         </td></tr>
</table>
</blockquote>

<p>
As a simple demonstration, see example <a href=\"modelica://Modelica.Blocks.Examples.Noise.UniformNoise\">Blocks.Examples.Noise.UniformNoise</a>.
In the next diagram, a simulation result is shown for samplePeriod=0.02 s and uniform distribution with
y_min=-1, y_max=3:
</p>
<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Examples/Noise/UniformNoise.png\">
</blockquote>

<h4>Advanced tab: General settings</h4>
<p>
In the <strong>Advanced</strong> tab of the parameter menu, further options can be set in the noise blocks
as shown in the next table:
</p>

<blockquote>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Parameter</th>
    <th>Description</th></tr>

<tr><td> enableNoise </td>
    <td> = true, if noise is generated at the output of the block (this is the default).<br>
         = false, if noise generation is switched off and the constant value
         y_off is provided as output.</td></tr>
<tr><td> y_off </td>
    <td> If enableNoise = false, the output of the block instance has the
         value y_off. Default is y_off = 0.0.
         Furthermore, if enableNoise = true and time&lt;startTime, the output of the block is also
         y_off (see description of parameter startTime below).</td></tr>
</table>
</blockquote>

<h4>Advanced tab: Initialization</h4>

<p>
For every block instance, the internally used pseudo random number generator
has its own state. This state must be properly initialized, depending on
the desired situation. For this purpose the following parameters can be defined:
</p>

<blockquote>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Parameter</th>
    <th>Description</th></tr>

<tr><td> useGlobalSeed </td>
    <td> = true, if the seed (= Integer number) defined in the \"inner GlobalSeed globalSeed\"
         component is used for the initialization of the random number generator used in this
         instance of the noise block.
         Therefore, whenever the globalSeed defines a different number, the noise at every
         instance is changing. This is the default setting and therefore the globalSeed component
         defines whether every new simulation run shall provide the same noise
         (e.g. for a parameter optimization of controller parameters), or
         whether every new simulation run shall provide different noise
         (e.g. for a Monte Carlo simulation).<br>
         = false, if the seed defined by globalSeed is ignored. For example, if
         aerodynamic turbulence is modelled with a noise block and this turbulence
         model shall be used for all simulation runs of a Monte Carlo simulation, then
         useGlobalSeed has to be set to false.</td></tr>

<tr><td> useAutomaticLocalSeed </td>
    <td> An Integer number, called local seed, is needed to initialize the random number
         generator for a specific block instance. Instances using the same local seed
         produce exactly the same random number values (so the same noise, if the other settings
         of the instances are the same).<br>
         If <strong>useAutomaticLocalSeed = true</strong>, the
         local seed is determined automatically using a hash value of the instance name of the model that is
         inquired with the Modelica built-in operator <a href=\"https://specification.modelica.org/v3.4/Ch3.html#getinstancename\">getInstanceName()</a>.
         Note, this means that the noise changes if the component is renamed.<br>
         If <strong>useAutomaticLocalSeed = false</strong>, the local seed is defined
         explicitly by parameter fixedLocalSeed. It is then guaranteed that the generated noise
         remains always the same (provided the other parameter values are the same).</td></tr>

<tr><td> fixedLocalSeed </td>
    <td> If useAutomaticLocalSeed = false, the local seed to be used.
         fixedLocalSeed can be any Integer number (including zero or a negative number).
         The initialization algorithm produces a meaningful initial state of the random
         number generator from fixedLocalSeed and (if useAutomaticGlobalSeed=true) from globalSeed even for
         bad seeds such as 0 or 1, so the subsequently drawing of random numbers produces always statistically
         meaningful numbers.</td></tr>

<tr><td> startTime </td>
    <td> The time instant at which noise shall be generated at the output y. The default
         startTime = 0.
         For time&lt;startTime, y = y_off. In some cases it is meaningful to simulate
         a certain duration until an approximate steady-state is reached. In such a case
         startTime should be set to a time instant after this duration.</td></tr>
</table>
</blockquote>

<h4>Random Number Generators</h4>

<p>
The core of the noise generation is the computation of uniform random
numbers in the range 0.0 .. 1.0 (and these random numbers are transformed
afterwards, see below). This sublibrary uses the xorshift random number generation
suite developed in 2014 by Sebastiano Vigna (for details see
<a href=\"http://xorshift.di.unimi.it\">http://xorshift.di.unimi.it</a> and
<a href=\"modelica://Modelica.Math.Random.Generators\">Math.Random.Generators</a>).
These random number generators have excellent
statistical properties, produce quickly statistically relevant random numbers, even if
starting from a bad initial seed, and have a reasonable length of the internal state
vector of 2, 4, and 33 Integer elements. The random number generator with an internal
state vector of length 2 is used to initialize the other two random number generators.
The length 4 random number generator is used in the noise blocks of this package, and every
such block has its own internal state vector, as needed for reproducible noise blocks.
The random number generator with a length of 33 Integer is used from the impure random number
generator. It is suited even for massively parallel simulations where every simulation
computes a large number of random values. More details of the random number
generators are described in the documentation of package
<a href=\"modelica://Modelica.Math.Random.Generators\">Math.Random.Generators</a>.
</p>

<h4>Distributions</h4>

<p>
The uniform random numbers in the range 0.0 .. 1.0 are transformed to a desired
random number distribution by selecting an appropriate <strong>distribution</strong> or
<strong>truncated distribution</strong>. For an example of a truncated distribution, see the following
diagram of the probability density function of a normal distribution
compared with its truncated version:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/TruncatedNormal.density.png\">
</blockquote>

<p>
The corresponding inverse cumulative distribution functions are shown in the next diagram:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/TruncatedNormal.quantile.png\">
</blockquote>

<p>
When providing an x-value between 0.0 .. 1.0 from a random number generator, then the truncated
inverse cumulative probability density function of a normal distribution transforms this value into the
desired band (in the diagram above to the range: -1.5 .. 1.5). Contrary to a standard distribution,
truncated distributions have the advantage that the resulting random values are guaranteed
to be in the defined band (whereas a standard normal distribution might also result in any value;
when modeling noise that is known to be in a particular range, say &plusmn; 0.1 Volt,
then with the TruncatedNormal distribution it is guaranteed that random values are only
generated in this band). More details of truncated
distributions are given in the documentation of package
<a href=\"modelica://Modelica.Math.Distributions\">Math.Distributions</a>.
</p>
</html>",     revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
    end Noise;

    package Sources
    "Library of signal source blocks generating Real, Integer and Boolean signals"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.SourcesPackage;

      block Constant "Generate constant signal of type Real"
        parameter Real k(start=1) "Constant output value"
        annotation(Dialog(groupImage="modelica://Modelica/Resources/Images/Blocks/Sources/Constant.png"));
        extends Interfaces.SO;

      equation
        y = k;
        annotation (
          defaultComponentName="const",
          Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-80,68},{-80,-80}}, color={192,192,192}),
              Polygon(
                points={{-80,90},{-88,68},{-72,68},{-80,90}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,-70},{82,-70}}, color={192,192,192}),
              Polygon(
                points={{90,-70},{68,-62},{68,-78},{90,-70}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,0},{80,0}}),
              Text(
                extent={{-150,-150},{150,-110}},
                textString="k=%k")}),
          Documentation(info="<html>
<p>
The Real output y is a constant signal:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/Constant.png\"
     alt=\"Constant.png\">
</p>
</html>"));
      end Constant;

      block BooleanConstant "Generate constant signal of type Boolean"
        parameter Boolean k=true "Constant output value"
        annotation(Dialog(groupImage="modelica://Modelica/Resources/Images/Blocks/Sources/BooleanConstant.png"));
        extends Interfaces.partialBooleanSource;

      equation
        y = k;
        annotation (
          Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={Line(points={{-80,0},{80,0}}),
                Text(
                extent={{-150,-140},{150,-110}},
                textString="%k")}),
            Documentation(info="<html>
<p>
The Boolean output y is a constant signal:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/BooleanConstant.png\"
     alt=\"BooleanConstant.png\">
</p>
</html>"));
      end BooleanConstant;
      annotation (Documentation(info="<html>
<p>
This package contains <strong>source</strong> components, i.e., blocks which
have only output signals. These blocks are used as signal generators
for Real, Integer and Boolean signals.
</p>

<p>
All Real source signals (with the exception of the Constant source)
have at least the following two parameters:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>offset</strong></td>
      <td>Value which is added to the signal</td>
  </tr>
  <tr><td><strong>startTime</strong></td>
      <td>Start time of signal. For time &lt; startTime,
                the output y is set to offset.</td>
  </tr>
</table>

<p>
The <strong>offset</strong> parameter is especially useful in order to shift
the corresponding source, such that at initial time the system
is stationary. To determine the corresponding value of offset,
usually requires a trimming calculation.
</p>
</html>",     revisions="<html>
<ul>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       Integer sources added. Step, TimeTable and BooleanStep slightly changed.</li>
<li><em>Nov. 8, 1999</em>
       by <a href=\"mailto:christoph@clauss-it.com\">Christoph Clau&szlig;</a>,
       <a href=\"mailto:Andre.Schneider@eas.iis.fraunhofer.de\">Andre.Schneider@eas.iis.fraunhofer.de</a>,
       <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       New sources: Exponentials, TimeTable. Trapezoid slightly enhanced
       (nperiod=-1 is an infinite number of periods).</li>
<li><em>Oct. 31, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       <a href=\"mailto:christoph@clauss-it.com\">Christoph Clau&szlig;</a>,
       <a href=\"mailto:Andre.Schneider@eas.iis.fraunhofer.de\">Andre.Schneider@eas.iis.fraunhofer.de</a>,
       All sources vectorized. New sources: ExpSine, Trapezoid,
       BooleanConstant, BooleanStep, BooleanPulse, SampleTrigger.
       Improved documentation, especially detailed description of
       signals in diagram layer.</li>
<li><em>June 29, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized a first version, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.</li>
</ul>
</html>"));
    end Sources;

    package Tables
    "Library of blocks to interpolate in one and two-dimensional tables"
      extends Modelica.Icons.Package;

      block CombiTable1Ds
        "Table look-up in one dimension (matrix/file) with one input and n outputs"
        extends Modelica.Blocks.Interfaces.SIMO(final nout=size(columns, 1));
        parameter Boolean tableOnFile=false
          "= true, if table is defined on file or in function usertab"
          annotation (Dialog(group="Table data definition"));
        parameter Real table[:, :] = fill(0.0, 0, 2)
          "Table matrix (grid = first column; e.g., table=[0, 0; 1, 1; 2, 4])"
          annotation (Dialog(group="Table data definition",enable=not tableOnFile));
        parameter String tableName="NoName"
          "Table name on file or in function usertab (see docu)"
          annotation (Dialog(group="Table data definition",enable=tableOnFile));
        parameter String fileName="NoName" "File where matrix is stored"
          annotation (Dialog(
            group="Table data definition",
            enable=tableOnFile,
            loadSelector(filter="Text files (*.txt);;MATLAB MAT-files (*.mat)",
                caption="Open file in which table is present")));
        parameter Boolean verboseRead=true
          "= true, if info message that file is loading is to be printed"
          annotation (Dialog(group="Table data definition",enable=tableOnFile));
        parameter Integer columns[:]=2:size(table, 2)
          "Columns of table to be interpolated"
          annotation (Dialog(group="Table data interpretation"));
        parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
          "Smoothness of table interpolation"
          annotation (Dialog(group="Table data interpretation"));
        parameter Modelica.Blocks.Types.Extrapolation extrapolation=Modelica.Blocks.Types.Extrapolation.LastTwoPoints
          "Extrapolation of data outside the definition range"
          annotation (Dialog(group="Table data interpretation"));
        parameter Boolean verboseExtrapolation=false
          "= true, if warning messages are to be printed if table input is outside the definition range"
          annotation (Dialog(group="Table data interpretation", enable=extrapolation == Modelica.Blocks.Types.Extrapolation.LastTwoPoints or extrapolation == Modelica.Blocks.Types.Extrapolation.HoldLastPoint));
        final parameter Real u_min=Internal.getTable1DAbscissaUmin(tableID)
          "Minimum abscissa value defined in table";
        final parameter Real u_max=Internal.getTable1DAbscissaUmax(tableID)
          "Maximum abscissa value defined in table";
      protected
        parameter Modelica.Blocks.Types.ExternalCombiTable1D tableID=
            Modelica.Blocks.Types.ExternalCombiTable1D(
              if tableOnFile then tableName else "NoName",
              if tableOnFile and fileName <> "NoName" and not Modelica.Utilities.Strings.isEmpty(fileName) then fileName else "NoName",
              table,
              columns,
              smoothness,
              extrapolation,
              if tableOnFile then verboseRead else false) "External table object";
      equation
        if tableOnFile then
          assert(tableName <> "NoName",
            "tableOnFile = true and no table name given");
        else
          assert(size(table, 1) > 0 and size(table, 2) > 0,
            "tableOnFile = false and parameter table is an empty matrix");
        end if;

        if verboseExtrapolation and (
          extrapolation == Modelica.Blocks.Types.Extrapolation.LastTwoPoints or
          extrapolation == Modelica.Blocks.Types.Extrapolation.HoldLastPoint) then
          assert(noEvent(u >= u_min), "
Extrapolation warning: The value u (="     + String(u) + ") must be greater or equal
than the minimum abscissa value u_min (="     + String(u_min) + ") defined in the table.
",     level=AssertionLevel.warning);
          assert(noEvent(u <= u_max), "
Extrapolation warning: The value u (="     + String(u) + ") must be less or equal
than the maximum abscissa value u_max (="     + String(u_max) + ") defined in the table.
",     level=AssertionLevel.warning);
        end if;

        if smoothness == Modelica.Blocks.Types.Smoothness.ConstantSegments then
          for i in 1:nout loop
            y[i] = Internal.getTable1DValueNoDer(tableID, i, u);
          end for;
        elseif smoothness == Modelica.Blocks.Types.Smoothness.LinearSegments then
          for i in 1:nout loop
            y[i] = Internal.getTable1DValueNoDer2(tableID, i, u);
          end for;
        else
          for i in 1:nout loop
            y[i] = Internal.getTable1DValue(tableID, i, u);
          end for;
        end if;
        annotation (
          Documentation(info="<html>
<p>
<strong>Univariate constant</strong>, <strong>linear</strong> or <strong>cubic Hermite
spline interpolation</strong> in <strong>one</strong> dimension of a
<strong>table</strong>.
Via parameter <strong>columns</strong> it can be defined how many columns of the
table are interpolated. If, e.g., columns={2,4}, it is assumed that
2 output signals are present and that the first output interpolates
via column 2 and the second output interpolates via column 4 of the
table matrix.
</p>
<p>
The grid points and function values are stored in a matrix \"table[i,j]\",
where the first column \"table[:,1]\" contains the grid points and the
other columns contain the data to be interpolated. Example:
</p>
<blockquote><pre>
table = [0,  0;
         1,  1;
         2,  4;
         4, 16]
If, e.g., the input u = 1.0, the output y =  1.0,
    e.g., the input u = 1.5, the output y =  2.5,
    e.g., the input u = 2.0, the output y =  4.0,
    e.g., the input u =-1.0, the output y = -1.0 (i.e., extrapolation).
</pre></blockquote>
<ul>
<li>The interpolation interval is found by a binary search where the interval used in the
    last call is used as start interval.</li>
<li>Via parameter <strong>smoothness</strong> it is defined how the data is interpolated:
<blockquote><pre>
smoothness = 1: Linear interpolation
           = 2: Akima interpolation: Smooth interpolation by cubic Hermite
                splines such that der(y) is continuous, also if extrapolated.
           = 3: Constant segments
           = 4: Fritsch-Butland interpolation: Smooth interpolation by cubic
                Hermite splines such that y preserves the monotonicity and
                der(y) is continuous, also if extrapolated.
           = 5: Steffen interpolation: Smooth interpolation by cubic Hermite
                splines such that y preserves the monotonicity and der(y)
                is continuous, also if extrapolated.
           = 6: Modified Akima interpolation: Smooth interpolation by cubic
                Hermite splines such that der(y) is continuous, also if
                extrapolated. Additionally, overshoots and edge cases of the
                original Akima interpolation method are avoided.
</pre></blockquote></li>
<li>First and second <strong>derivatives</strong> are provided, with exception of the following two smoothness options.
<ol>
<li>No derivatives are provided for interpolation by constant segments.</li>
<li>No second derivative is provided for linear interpolation.</li>
</ol></li>
<li>Values <strong>outside</strong> of the table range, are computed by
    extrapolation according to the setting of parameter <strong>extrapolation</strong>:
<blockquote><pre>
extrapolation = 1: Hold the first or last value of the table,
                   if outside of the table scope.
              = 2: Extrapolate by using the derivative at the first/last table
                   points if outside of the table scope.
                   (If smoothness is LinearSegments or ConstantSegments
                   this means to extrapolate linearly through the first/last
                   two table points.).
              = 3: Periodically repeat the table data (periodical function).
              = 4: No extrapolation, i.e. extrapolation triggers an error
</pre></blockquote></li>
<li>If the table has only <strong>one row</strong>, the table value is returned,
    independent of the value of the input signal.</li>
<li>The grid values (first column) have to be strictly increasing.</li>
</ul>
<p>
The table matrix can be defined in the following ways:
</p>
<ol>
<li>Explicitly supplied as <strong>parameter matrix</strong> \"table\",
    and the other parameters have the following values:
<blockquote><pre>
tableName is \"NoName\" or has only blanks,
fileName  is \"NoName\" or has only blanks.
</pre></blockquote></li>
<li><strong>Read</strong> from a <strong>file</strong> \"fileName\" where the matrix is stored as
    \"tableName\". Both text and MATLAB MAT-file format is possible.
    (The text format is described below).
    The MAT-file format comes in four different versions: v4, v6, v7 and v7.3.
    The library supports at least v4, v6 and v7 whereas v7.3 is optional.
    It is most convenient to generate the MAT-file from FreeMat or MATLAB&reg;
    by command
<blockquote><pre>
save tables.mat tab1 tab2 tab3
</pre></blockquote>
    or Scilab by command
<blockquote><pre>
savematfile tables.mat tab1 tab2 tab3
</pre></blockquote>
    when the three tables tab1, tab2, tab3 should be used from the model.<br>
    Note, a fileName can be defined as URI by using the helper function
    <a href=\"modelica://Modelica.Utilities.Files.loadResource\">loadResource</a>.</li>
<li>Statically stored in function \"usertab\" in file \"usertab.c\".
    The matrix is identified by \"tableName\". Parameter
    fileName = \"NoName\" or has only blanks. Row-wise storage is always to be
    preferred as otherwise the table is reallocated and transposed.
    See the <a href=\"modelica://Modelica.Blocks.Tables\">Tables</a> package
    documentation for more details.</li>
</ol>
<p>
When the constant \"NO_FILE_SYSTEM\" is defined, all file I/O related parts of the
source code are removed by the C-preprocessor, such that no access to files takes place.
</p>
<p>
If tables are read from a text file, the file needs to have the
following structure (\"-----\" is not part of the file content):
</p>
<blockquote><pre>
-----------------------------------------------------
#1
double tab1(5,2)   # comment line
  0   0
  1   1
  2   4
  3   9
  4  16
double tab2(5,2)   # another comment line
  0   0
  2   2
  4   8
  6  18
  8  32
-----------------------------------------------------
</pre></blockquote>
<p>
Note, that the first two characters in the file need to be
\"#1\" (a line comment defining the version number of the file format).
Afterwards, the corresponding matrix has to be declared
with type (= \"double\" or \"float\"), name and actual dimensions.
Finally, in successive rows of the file, the elements of the matrix
have to be given. The elements have to be provided as a sequence of
numbers in row-wise order (therefore a matrix row can span several
lines in the file and need not start at the beginning of a line).
Numbers have to be given according to C syntax (such as 2.3, -2, +2.e4).
Number separators are spaces, tab (\\t), comma (,), or semicolon (;).
Several matrices may be defined one after another. Line comments start
with the hash symbol (#) and can appear everywhere.
Text files should either be ASCII or UTF-8 encoded, where UTF-8 encoded strings are only allowed in line comments and an optional UTF-8 BOM at the start of the text file is ignored.
Other characters, like trailing non comments, are not allowed in the file.
</p>
<p>
MATLAB is a registered trademark of The MathWorks, Inc.
</p>
</html>"),Icon(
          coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Line(points={{-60.0,40.0},{-60.0,-40.0},{60.0,-40.0},{60.0,40.0},{30.0,40.0},{30.0,-40.0},{-30.0,-40.0},{-30.0,40.0},{-60.0,40.0},{-60.0,20.0},{60.0,20.0},{60.0,0.0},{-60.0,0.0},{-60.0,-20.0},{60.0,-20.0},{60.0,-40.0},{-60.0,-40.0},{-60.0,40.0},{60.0,40.0},{60.0,-40.0}}),
          Line(points={{0.0,40.0},{0.0,-40.0}}),
          Rectangle(fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-60.0,20.0},{-30.0,40.0}}),
          Rectangle(fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-60.0,0.0},{-30.0,20.0}}),
          Rectangle(fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-60.0,-20.0},{-30.0,0.0}}),
          Rectangle(fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-60.0,-40.0},{-30.0,-20.0}})}));
      end CombiTable1Ds;

      package Internal
      "Internal external object definitions for table functions that should not be directly utilized by the user"
        extends Modelica.Icons.InternalPackage;

        pure function getTable1DValue "Interpolate 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          input Integer icol "Column number";
          input Real u "Abscissa value";
          output Real y "Interpolated value";
          external "C" y = ModelicaStandardTables_CombiTable1D_getValue(tableID, icol, u)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
          annotation (derivative = getDerTable1DValue);
        end getTable1DValue;

        pure function getTable1DValueNoDer
          "Interpolate 1-dim. table defined by matrix (but do not provide a derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          input Integer icol "Column number";
          input Real u "Abscissa value";
          output Real y "Interpolated value";
          external "C" y = ModelicaStandardTables_CombiTable1D_getValue(tableID, icol, u)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTable1DValueNoDer;

        pure function getTable1DValueNoDer2
          "Interpolate 1-dim. table defined by matrix (but do not provide a second derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          input Integer icol "Column number";
          input Real u "Abscissa value";
          output Real y "Interpolated value";
          external "C" y = ModelicaStandardTables_CombiTable1D_getValue(tableID, icol, u)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
          annotation (derivative = getDerTable1DValueNoDer);
        end getTable1DValueNoDer2;

        pure function getDerTable1DValue
          "Derivative of interpolated 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          input Integer icol "Column number";
          input Real u "Abscissa value";
          input Real der_u "Derivative of abscissa value";
          output Real der_y "Derivative of interpolated value";
          external "C" der_y = ModelicaStandardTables_CombiTable1D_getDerValue(tableID, icol, u, der_u)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
          annotation (derivative(order=2) = getDer2Table1DValue);
        end getDerTable1DValue;

        pure function getDerTable1DValueNoDer
          "Derivative of interpolated 1-dim. table defined by matrix (but do not provide a second derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          input Integer icol "Column number";
          input Real u "Abscissa value";
          input Real der_u "Derivative of abscissa value";
          output Real der_y "Derivative of interpolated value";
          external "C" der_y = ModelicaStandardTables_CombiTable1D_getDerValue(tableID, icol, u, der_u)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getDerTable1DValueNoDer;

        pure function getDer2Table1DValue
          "Second derivative of interpolated 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          input Integer icol "Column number";
          input Real u "Abscissa value";
          input Real der_u "Derivative of abscissa value";
          input Real der2_u " Second derivative of abscissa value";
          output Real der2_y "Second derivative of interpolated value";
          external "C" der2_y = ModelicaStandardTables_CombiTable1D_getDer2Value(tableID, icol, u, der_u, der2_u)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getDer2Table1DValue;

        pure function getTable1DAbscissaUmin
          "Return minimum abscissa value of 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          output Real uMin "Minimum abscissa value in table";
          external "C" uMin = ModelicaStandardTables_CombiTable1D_minimumAbscissa(tableID)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTable1DAbscissaUmin;

        pure function getTable1DAbscissaUmax
          "Return maximum abscissa value of 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable1D tableID "External table object";
          output Real uMax "Maximum abscissa value in table";
          external "C" uMax = ModelicaStandardTables_CombiTable1D_maximumAbscissa(tableID)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTable1DAbscissaUmax;
      end Internal;
      annotation (Documentation(info="<html>
<p>This package contains blocks for one- and two-dimensional interpolation in tables.</p>
<h4>Special interest topic: Statically stored tables for real-time simulation targets</h4>
<p>Especially for use on real-time platform targets (e.g., HIL-simulators) with <strong>no file system</strong>, it is possible to statically
store tables using a function &quot;usertab&quot; in a file conventionally named &quot;usertab.c&quot;. This can be more efficient than providing the tables as Modelica parameter arrays.</p>
<p>This is achieved by providing the tables in a specific structure as C-code and compiling that C-code together with the rest of the simulation model into a binary
that can be executed on the target platform. The &quot;Resources/Data/Tables/&quot; subdirectory of the MSL installation directory contains the files
<a href=\"modelica://Modelica/Resources/Data/Tables/usertab.c\">&quot;usertab.c&quot;</a> and <a href=\"modelica://Modelica/Resources/Data/Tables/usertab.h\">&quot;usertab.h&quot;</a>
that can be used as a template for own developments. While &quot;usertab.c&quot; would be typically used unmodified, the
&quot;usertab.h&quot; needs to adapted for the own needs.</p>
<p>In order to work it is necessary that the compiler pulls in the &quot;usertab.c&quot; file. Different Modelica tools might provide different mechanisms to do so.
Please consult the respective documentation/support for your Modelica tool.</p>
<p>A possible (though slightly makeshift) approach is to pull in the required files by utilizing a &quot;dummy&quot;-function that uses the Modelica external function
interface to include the required &quot;usertab.c&quot;. An example how this can be done is given below.</p>
<blockquote><pre>
model ExampleCTable \"Example utilizing the usertab.c interface\"
  extends Modelica.Icons.Example;
  parameter Real dummy(fixed=false) \"Dummy parameter\" annotation(HideResult=true);
  Modelica.Blocks.Tables.CombiTable1Dv table(tableOnFile=true, tableName=\"TestTable_1D_a\")
    annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Modelica.Blocks.Sources.ContinuousClock clock
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
protected
  encapsulated impure function getUsertab \"External dummy function to include \\\"usertab.c\\\"\"
    input Real dummy_u[:];
    output Real dummy_y;
    external \"C\" dummy_y = mydummyfunc(dummy_u);
    annotation(IncludeDirectory=\"modelica://Modelica/Resources/Data/Tables\",
           Include = \"#include \"usertab.c\"
double mydummyfunc(double* dummy_in) {
   return 0;
}
\");
  end getUsertab;
initial equation
  dummy = getUsertab(table.y);
equation
  connect(clock.y, table.u[1]) annotation (Line(points={{-59,10},{-42,10}}, color={0,0,127}));
  annotation (experiment(StartTime=0, StopTime=5), uses(Modelica(version=\"4.0.0\")));
end ExampleCTable;
</pre></blockquote>
</html>"),     Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Rectangle(
              extent={{-76,-26},{80,-76}},
              lineColor={95,95,95},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-76,24},{80,-26}},
              lineColor={95,95,95},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-76,74},{80,24}},
              lineColor={95,95,95},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-28,74},{-28,-76}},
              color={95,95,95}),
            Line(
              points={{24,74},{24,-76}},
              color={95,95,95})}));
    end Tables;

    package Types
    "Library of constants, external objects and types with choices, especially to build menus"
      extends Modelica.Icons.TypesPackage;

      type Smoothness = enumeration(
          LinearSegments "Linear interpolation of table points",
          ContinuousDerivative
            "Akima spline interpolation of table points (such that the first derivative is continuous)",
          ConstantSegments
            "Piecewise constant interpolation of table points (the value from the previous abscissa point is returned)",
          MonotoneContinuousDerivative1
            "Fritsch-Butland spline interpolation (such that the monotonicity is preserved and the first derivative is continuous)",
          MonotoneContinuousDerivative2
            "Steffen spline interpolation of table points (such that the monotonicity is preserved and the first derivative is continuous)",
          ModifiedContinuousDerivative
            "Modified Akima spline interpolation of table points (such that the first derivative is continuous and shortcomings of the original Akima method are avoided)")
        "Enumeration defining the smoothness of table interpolation";

        type Extrapolation = enumeration(
          HoldLastPoint
            "Hold the first/last table point outside of the table scope",
          LastTwoPoints
            "Extrapolate by using the derivative at the first/last table points outside of the table scope",
          Periodic "Repeat the table scope periodically",
          NoExtrapolation "Extrapolation triggers an error")
        "Enumeration defining the extrapolation of table interpolation";

      class ExternalCombiTable1D
        "External object of 1-dim. table defined by matrix"
        extends ExternalObject;

        function constructor "Initialize 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input String tableName "Table name";
          input String fileName "File name";
          input Real table[:, :];
          input Integer columns[:];
          input Modelica.Blocks.Types.Smoothness smoothness;
          input Modelica.Blocks.Types.Extrapolation extrapolation=Modelica.Blocks.Types.Extrapolation.LastTwoPoints;
          input Boolean verboseRead=true "= true: Print info message; = false: No info message";
          output ExternalCombiTable1D externalCombiTable1D;
        external "C" externalCombiTable1D = ModelicaStandardTables_CombiTable1D_init2(
                fileName,
                tableName,
                table,
                size(table, 1),
                size(table, 2),
                columns,
                size(columns, 1),
                smoothness,
                extrapolation,
                verboseRead) annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end constructor;

        function destructor "Terminate 1-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input ExternalCombiTable1D externalCombiTable1D;
        external "C" ModelicaStandardTables_CombiTable1D_close(externalCombiTable1D)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end destructor;

      end ExternalCombiTable1D;
      annotation (Documentation(info="<html>
<p>
In this package <strong>types</strong>, <strong>constants</strong> and <strong>external objects</strong> are defined that are used
in library Modelica.Blocks. The types have additional annotation choices
definitions that define the menus to be built up in the graphical
user interface when the type is used as parameter in a declaration.
</p>
</html>"));
    end Types;

    package Icons "Icons for Blocks"
        extends Modelica.Icons.IconsPackage;

        partial block Block "Basic graphical layout of input/output block"

          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={Rectangle(
                extent={{-100,-100},{100,100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-150,150},{150,110}},
                textString="%name",
                textColor={0,0,255})}),
          Documentation(info="<html>
<p>
Block that has only the basic icon for an input/output
block (no declarations, no equations). Most blocks
of package Modelica.Blocks inherit directly or indirectly
from this block.
</p>
</html>"));

        end Block;

      partial block PartialBooleanBlock "Basic graphical layout of logical block"

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={210,210,210},
                fillPattern=FillPattern.Solid,
                borderPattern=BorderPattern.Raised), Text(
                extent={{-150,150},{150,110}},
                textString="%name",
                textColor={0,0,255})}), Documentation(info="<html>
<p>
Block that has only the basic icon for an input/output,
Boolean block (no declarations, no equations) used especially
in the Blocks.Logical library.
</p>
</html>"));
      end PartialBooleanBlock;
    end Icons;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
        Rectangle(
          origin={0.0,35.1488},
          fillColor={255,255,255},
          extent={{-30.0,-20.1488},{30.0,20.1488}}),
        Rectangle(
          origin={0.0,-34.8512},
          fillColor={255,255,255},
          extent={{-30.0,-20.1488},{30.0,20.1488}}),
        Line(
          origin={-51.25,0.0},
          points={{21.25,-35.0},{-13.75,-35.0},{-13.75,35.0},{6.25,35.0}}),
        Polygon(
          origin={-40.0,35.0},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{10.0,0.0},{-5.0,5.0},{-5.0,-5.0}}),
        Line(
          origin={51.25,0.0},
          points={{-21.25,35.0},{13.75,35.0},{13.75,-35.0},{-6.25,-35.0}}),
        Polygon(
          origin={40.0,-35.0},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{-10.0,0.0},{5.0,5.0},{5.0,-5.0}})}), Documentation(info="<html>
<p>
This library contains input/output blocks to build up block diagrams.
</p>

<dl>
<dt><strong>Main Author:</strong></dt>
<dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
    Deutsches Zentrum f&uuml;r Luft und Raumfahrt e. V. (DLR)<br>
    Oberpfaffenhofen<br>
    Postfach 1116<br>
    D-82230 Wessling<br>
    email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a><br></dd>
</dl>
<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>June 23, 2004</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Introduced new block connectors and adapted all blocks to the new connectors.
       Included subpackages Continuous, Discrete, Logical, Nonlinear from
       package ModelicaAdditions.Blocks.
       Included subpackage ModelicaAdditions.Table in Modelica.Blocks.Sources
       and in the new package Modelica.Blocks.Tables.
       Added new blocks to Blocks.Sources and Blocks.Logical.
       </li>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       New subpackage Examples, additional components.
       </li>
<li><em>June 20, 2000</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and
       Michael Tiller:<br>
       Introduced a replaceable signal type into
       Blocks.Interfaces.RealInput/RealOutput:
<blockquote><pre>
replaceable type SignalType = Real
</pre></blockquote>
       in order that the type of the signal of an input/output block
       can be changed to a physical type, for example:
<blockquote><pre>
Sine sin1(outPort(redeclare type SignalType=Modelica.Units.SI.Torque))
</pre></blockquote>
      </li>
<li><em>Sept. 18, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Renamed to Blocks. New subpackages Math, Nonlinear.
       Additional components in subpackages Interfaces, Continuous
       and Sources.</li>
<li><em>June 30, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized a first version, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.</li>
</ul>
</html>"));
  end Blocks;

  package Electrical
  "Library of electrical models (analog, digital, machines, polyphase)"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;

    package Analog "Library for analog electrical models"
      extends Modelica.Icons.Package;

      package Interfaces
      "Connectors and partial models for Analog electrical components"
        extends Modelica.Icons.InterfacesPackage;

        connector Pin "Pin of an electrical component"
          SI.ElectricPotential v "Potential at the pin" annotation (
              unassignedMessage="An electrical potential cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
          flow SI.Current i "Current flowing into the pin" annotation (
              unassignedMessage="An electrical current cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
          annotation (defaultComponentName="pin",
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                    100}}), graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={Rectangle(
                  extent={{-40,40},{40,-40}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-160,110},{40,50}},
                  textColor={0,0,255},
                  textString="%name")}),
            Documentation(revisions="<html>
<ul>
<li><em> 1998   </em>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>",       info="<html>
<p>Pin is the basic electric connector. It includes the voltage which consists between the pin and the ground node. The ground node is the node of (any) ground device (Modelica.Electrical.Basic.Ground). Furthermore, the pin includes the current, which is considered to be <strong>positive</strong> if it is flowing at the pin <strong>into the device</strong>.</p>
</html>"));
        end Pin;

        connector PositivePin "Positive pin of an electrical component"
          SI.ElectricPotential v "Potential at the pin" annotation (
              unassignedMessage="An electrical potential cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
          flow SI.Current i "Current flowing into the pin" annotation (
              unassignedMessage="An electrical current cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
          annotation (defaultComponentName="pin_p",
            Documentation(info="<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>",                     revisions="<html>
<ul>
<li><em> 1998   </em>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>"),  Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                    100}}), graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={Rectangle(
                  extent={{-40,40},{40,-40}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-160,110},{40,50}},
                  textColor={0,0,255},
                  textString="%name")}));
        end PositivePin;

        connector NegativePin "Negative pin of an electrical component"
          SI.ElectricPotential v "Potential at the pin" annotation (
              unassignedMessage="An electrical potential cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
          flow SI.Current i "Current flowing into the pin" annotation (
              unassignedMessage="An electrical current cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
          annotation (defaultComponentName="pin_n",
            Documentation(info="<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>",       revisions="<html>
<dl>
<dt><em>1998</em></dt>
<dd>by Christoph Clauss initially implemented
</dd>
</dl>
</html>"),  Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                    100}}), graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={Rectangle(
                  extent={{-40,40},{40,-40}},
                  lineColor={0,0,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-40,110},{160,50}},
                  textString="%name",
                  textColor={0,0,255})}));
        end NegativePin;
        annotation (Documentation(info="<html>
<p>This package contains connectors and interfaces (partial models) for analog electrical components. The partial models contain typical combinations of pins, and internal variables which are often used. Furthermore, the thermal heat port is in this package which can be included by inheritance.</p>
</html>",revisions="<html>
<dl>
<dt>
<strong>Main Authors:</strong>
</dt>
<dd>
Christoph Clau&szlig;
    &lt;<a href=\"mailto:christoph@clauss-it.com\">christoph@clauss-it.com</a>&gt;<br>
    Andr&eacute; Schneider
    &lt;<a href=\"mailto:Andre.Schneider@eas.iis.fraunhofer.de\">Andre.Schneider@eas.iis.fraunhofer.de</a>&gt;<br>
    Fraunhofer Institute for Integrated Circuits<br>
    Design Automation Department<br>
    Zeunerstra&szlig;e 38<br>
    D-01069 Dresden
</dd>
</dl>

<ul>
<li><em> 1998</em>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>"));
      end Interfaces;
      annotation (Documentation(info="<html>
<p>
This package contains packages for single-phase electrical components, see
<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide\">User&#39;s Guide</a></p>
</html>"),     Icon(graphics={
            Line(
              points={{12,60},{12,-60}}),
            Line(
              points={{-12,60},{-12,-60}}),
            Line(points={{-80,0},{-12,0}}),
            Line(points={{12,0},{80,0}})}));
    end Analog;
    annotation (
    Documentation(info="<html>
<p>
This library contains electrical components to build up analog and digital circuits,
as well as machines to model electrical motors and generators,
especially three-phase induction machines such as an asynchronous motor.
</p>

</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
      Rectangle(
        origin={20.3125,82.8571},
        extent={{-45.3125,-57.8571},{4.6875,-27.8571}}),
      Line(
        origin={8.0,48.0},
        points={{32.0,-58.0},{72.0,-58.0}}),
      Line(
        origin={9.0,54.0},
        points={{31.0,-49.0},{71.0,-49.0}}),
      Line(
        origin={-2.0,55.0},
        points={{-83.0,-50.0},{-33.0,-50.0}}),
      Line(
        origin={-3.0,45.0},
        points={{-72.0,-55.0},{-42.0,-55.0}}),
      Line(
        origin={1.0,50.0},
        points={{-61.0,-45.0},{-61.0,-10.0},{-26.0,-10.0}}),
      Line(
        origin={7.0,50.0},
        points={{18.0,-10.0},{53.0,-10.0},{53.0,-45.0}}),
      Line(
        origin={6.2593,48.0},
        points={{53.7407,-58.0},{53.7407,-93.0},{-66.2593,-93.0},{-66.2593,-58.0}})}));
  end Electrical;

  package Fluid
  "Library of 1-dim. thermo-fluid flow models using the Modelica.Media media description"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;
    import Cv = Modelica.Units.Conversions;

    model System
      "System properties and default values (ambient, flow direction, initialization)"

      package Medium = Modelica.Media.Interfaces.PartialMedium
        "Medium model for default start values"
          annotation (choicesAllMatching = true);
      parameter SI.AbsolutePressure p_ambient=101325
        "Default ambient pressure"
        annotation(Dialog(group="Environment"));
      parameter SI.Temperature T_ambient=293.15
        "Default ambient temperature"
        annotation(Dialog(group="Environment"));
      parameter SI.Acceleration g=Modelica.Constants.g_n
        "Constant gravity acceleration"
        annotation(Dialog(group="Environment"));

      // Assumptions
      parameter Boolean allowFlowReversal = true
        "= false to restrict to design flow direction (port_a -> port_b)"
        annotation(Dialog(tab="Assumptions"), Evaluate=true);
      parameter Modelica.Fluid.Types.Dynamics energyDynamics=
        Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
        "Default formulation of energy balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      parameter Modelica.Fluid.Types.Dynamics massDynamics=
        energyDynamics "Default formulation of mass balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      final parameter Modelica.Fluid.Types.Dynamics substanceDynamics=
        massDynamics "Default formulation of substance balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      final parameter Modelica.Fluid.Types.Dynamics traceDynamics=
        massDynamics "Default formulation of trace substance balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      parameter Modelica.Fluid.Types.Dynamics momentumDynamics=
        Modelica.Fluid.Types.Dynamics.SteadyState
        "Default formulation of momentum balances, if options available"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

      // Initialization
      parameter SI.MassFlowRate m_flow_start = 0
        "Default start value for mass flow rates"
        annotation(Dialog(tab = "Initialization"));
      parameter SI.AbsolutePressure p_start = p_ambient
        "Default start value for pressures"
        annotation(Dialog(tab = "Initialization"));
      parameter SI.Temperature T_start = T_ambient
        "Default start value for temperatures"
        annotation(Dialog(tab = "Initialization"));
      // Advanced
      parameter Boolean use_eps_Re = false
        "= true to determine turbulent region automatically using Reynolds number"
        annotation(Evaluate=true, Dialog(tab = "Advanced"));
      parameter SI.MassFlowRate m_flow_nominal = if use_eps_Re then 1 else 1e2*m_flow_small
        "Default nominal mass flow rate"
        annotation(Dialog(tab="Advanced", enable = use_eps_Re));
      parameter Real eps_m_flow(min=0) = 1e-4
        "Regularization of zero flow for |m_flow| < eps_m_flow*m_flow_nominal"
        annotation(Dialog(tab = "Advanced", enable = use_eps_Re));
      parameter SI.AbsolutePressure dp_small(min=0) = 1
        "Default small pressure drop for regularization of laminar and zero flow"
        annotation(Dialog(tab="Advanced", group="Classic", enable = not use_eps_Re));
      parameter SI.MassFlowRate m_flow_small(min=0) = 1e-2
        "Default small mass flow rate for regularization of laminar and zero flow"
        annotation(Dialog(tab = "Advanced", group="Classic", enable = not use_eps_Re));
    initial equation
      //assert(use_eps_Re, "*** Using classic system.m_flow_small and system.dp_small."
      //       + " They do not distinguish between laminar flow and regularization of zero flow."
      //       + " Absolute small values are error prone for models with local nominal values."
      //       + " Moreover dp_small can generally be obtained automatically."
      //       + " Please update the model to new system.use_eps_Re = true  (see system, Advanced tab). ***",
      //       level=AssertionLevel.warning);
      annotation (
        defaultComponentName="system",
        defaultComponentPrefixes="inner",
        missingInnerMessage="
Your model is using an outer \"system\" component but
an inner \"system\" component is not defined.
For simulation drag Modelica.Fluid.System into your model
to specify system properties.",
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textColor={0,0,255},
              textString="%name"),
            Line(points={{-86,-30},{82,-30}}),
            Line(points={{-82,-68},{-52,-30}}),
            Line(points={{-48,-68},{-18,-30}}),
            Line(points={{-14,-68},{16,-30}}),
            Line(points={{22,-68},{52,-30}}),
            Line(points={{74,84},{74,14}}),
            Polygon(
              points={{60,14},{88,14},{74,-18},{60,14}},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{16,20},{60,-18}},
              textString="g"),
            Text(
              extent={{-90,82},{74,50}},
              textString="defaults"),
            Line(
              points={{-82,14},{-42,-20},{2,30}},
              thickness=0.5),
            Ellipse(
              extent={{-10,40},{12,18}},
              pattern=LinePattern.None,
              fillColor={255,0,0},
              fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p>
 A system component is needed in each fluid model to provide system-wide settings, such as ambient conditions and overall modeling assumptions.
 The system settings are propagated to the fluid models using the inner/outer mechanism.
</p>
<p>
 A model should never directly use system parameters.
 Instead a local parameter should be declared, which uses the global setting as default.
 The only exceptions are:</p>
 <ul>
  <li>the gravity system.g,</li>
  <li>the global system.eps_m_flow, which is used to define a local m_flow_small for the local m_flow_nominal:
      <blockquote><pre>m_flow_small = system.eps_m_flow*m_flow_nominal</pre></blockquote>
  </li>
 </ul>
<p>
 The global system.m_flow_small and system.dp_small are classic parameters.
 They do not distinguish between laminar flow and regularization of zero flow.
 Absolute small values are error prone for models with local nominal values.
 Moreover dp_small can generally be obtained automatically.
 Consider using the new system.use_eps_Re = true (see Advanced tab).
</p>
</html>"));
    end System;

    package Interfaces
    "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow"
      extends Modelica.Icons.InterfacesPackage;

      connector FluidPort
        "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"

        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
          "Medium model" annotation (choicesAllMatching=true);

        flow Medium.MassFlowRate m_flow
          "Mass flow rate from the connection point into the component";
        Medium.AbsolutePressure p "Thermodynamic pressure in the connection point";
        stream Medium.SpecificEnthalpy h_outflow
          "Specific thermodynamic enthalpy close to the connection point if m_flow < 0";
        stream Medium.MassFraction Xi_outflow[Medium.nXi]
          "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
        stream Medium.ExtraProperty C_outflow[Medium.nC]
          "Properties c_i/m close to the connection point if m_flow < 0";
      end FluidPort;

      connector FluidPort_a "Generic fluid connector at design inlet"
        extends FluidPort;
        annotation (defaultComponentName="port_a",
                    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Ellipse(
                extent={{-40,40},{40,-40}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid), Text(extent={{-150,110},{150,50}},
                  textString="%name")}),
             Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,127,255},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid), Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end FluidPort_a;
      annotation (Documentation(info="<html>

</html>",     revisions="<html>
<ul>
<li><em>June 9th, 2008</em>
       by Michael Sielemann: Introduced stream keyword after decision at 57th Design Meeting (Lund).</li>
<li><em>May 30, 2007</em>
       by Christoph Richter: moved everything back to its original position in Modelica.Fluid.</li>
<li><em>Apr. 20, 2007</em>
       by Christoph Richter: moved parts of the original package from Modelica.Fluid
       to the development branch of Modelica 2.2.2.</li>
<li><em>Nov. 2, 2005</em>
       by Francesco Casella: restructured after 45th Design Meeting.</li>
<li><em>Nov. 20-21, 2002</em>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.</li>
<li><em>Nov. 11, 2002</em>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><em>Nov. 6, 2002</em>
       by Hilding Elmqvist: first version.</li>
<li><em>Aug. 11, 2002</em>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</ul>
</html>"));
    end Interfaces;

    package Types "Common types for fluid models"
      extends Modelica.Icons.TypesPackage;

      type Dynamics = enumeration(
          DynamicFreeInitial
            "DynamicFreeInitial -- Dynamic balance, Initial guess value",
          FixedInitial "FixedInitial -- Dynamic balance, Initial value fixed",
          SteadyStateInitial
            "SteadyStateInitial -- Dynamic balance, Steady state initial with guess value",
          SteadyState "SteadyState -- Steady state balance, Initial guess value")
        "Enumeration to define definition of balance equations"
      annotation (Documentation(info="<html>
<p>
Enumeration to define the formulation of balance equations
(to be selected via choices menu):
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th><strong>Dynamics.</strong></th><th><strong>Meaning</strong></th></tr>
<tr><td>DynamicFreeInitial</td><td>Dynamic balance, Initial guess value</td></tr>

<tr><td>FixedInitial</td><td>Dynamic balance, Initial value fixed</td></tr>

<tr><td>SteadyStateInitial</td><td>Dynamic balance, Steady state initial with guess value</td></tr>

<tr><td>SteadyState</td><td>Steady state balance, Initial guess value</td></tr>
</table>

<p>
The enumeration \"Dynamics\" is used for the mass, energy and momentum balance equations
respectively. The exact meaning for the three balance equations is stated in the following
tables:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><td colspan=\"3\"><strong>Mass balance</strong> </td></tr>
<tr><td><strong>Dynamics.</strong></td>
    <td><strong>Balance equation</strong></td>
    <td><strong>Initial condition</strong></td></tr>

<tr><td> DynamicFreeInitial</td>
    <td> no restrictions </td>
    <td> no initial conditions </td></tr>

<tr><td> FixedInitial</td>
    <td> no restrictions </td>
    <td> <strong>if</strong> Medium.singleState <strong>then</strong><br>
         &nbsp;&nbsp;no initial condition<br>
         <strong>else</strong> p=p_start </td></tr>

<tr><td> SteadyStateInitial</td>
    <td> no restrictions </td>
    <td> <strong>if</strong> Medium.singleState <strong>then</strong><br>
         &nbsp;&nbsp;no initial condition<br>
         <strong>else</strong> <strong>der</strong>(p)=0 </td></tr>

<tr><td> SteadyState</td>
    <td> <strong>der</strong>(m)=0  </td>
    <td> no initial conditions </td></tr>
</table>

&nbsp;<br>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><td colspan=\"3\"><strong>Energy balance</strong> </td></tr>
<tr><td><strong>Dynamics.</strong></td>
    <td><strong>Balance equation</strong></td>
    <td><strong>Initial condition</strong></td></tr>

<tr><td> DynamicFreeInitial</td>
    <td> no restrictions </td>
    <td> no initial conditions </td></tr>

<tr><td> FixedInitial</td>
    <td> no restrictions </td>
    <td> T=T_start or h=h_start </td></tr>

<tr><td> SteadyStateInitial</td>
    <td> no restrictions </td>
    <td> <strong>der</strong>(T)=0 or <strong>der</strong>(h)=0 </td></tr>

<tr><td> SteadyState</td>
    <td> <strong>der</strong>(U)=0  </td>
    <td> no initial conditions </td></tr>
</table>

&nbsp;<br>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><td colspan=\"3\"><strong>Momentum balance</strong> </td></tr>
<tr><td><strong>Dynamics.</strong></td>
    <td><strong>Balance equation</strong></td>
    <td><strong>Initial condition</strong></td></tr>

<tr><td> DynamicFreeInitial</td>
    <td> no restrictions </td>
    <td> no initial conditions </td></tr>

<tr><td> FixedInitial</td>
    <td> no restrictions </td>
    <td> m_flow = m_flow_start </td></tr>

<tr><td> SteadyStateInitial</td>
    <td> no restrictions </td>
    <td> <strong>der</strong>(m_flow)=0 </td></tr>

<tr><td> SteadyState</td>
    <td> <strong>der</strong>(m_flow)=0 </td>
    <td> no initial conditions </td></tr>
</table>

<p>
In the tables above, the equations are given for one-substance fluids. For multiple-substance
fluids and for trace substances, equivalent equations hold.
</p>

<p>
Medium.singleState is a medium property and defines whether the medium is only
described by one state (+ the mass fractions in case of a multi-substance fluid). In such
a case one initial condition less must be provided. For example, incompressible
media have Medium.singleState = <strong>true</strong>.
</p>

</html>"));
      annotation (preferredView="info",
                  Documentation(info="<html>

</html>"));
    end Types;
  annotation (Icon(graphics={
          Polygon(points={{-70,26},{68,-44},{68,26},{2,-10},{-70,-42},{-70,26}}),
          Line(points={{2,42},{2,-10}}),
          Rectangle(
            extent={{-18,50},{22,42}},
            fillPattern=FillPattern.Solid)}), preferredView="info",
    Documentation(info="<html>
<p>
Library <strong>Modelica.Fluid</strong> is a <strong>free</strong> Modelica package providing components for
<strong>1-dimensional thermo-fluid flow</strong> in networks of vessels, pipes, fluid machines, valves and fittings.
A unique feature is that the component equations and the media models
as well as pressure loss and heat transfer correlations are decoupled from each other.
All components are implemented such that they can be used for
media from the Modelica.Media library. This means especially that an
incompressible or compressible medium, a single or a multiple
substance medium with one or more phases might be used.
</p>

<p>
In the next figure, several features of the library are demonstrated with
a simple heating system with a closed flow cycle. By just changing one configuration parameter in the system object the equations are changed between steady-state and dynamic simulation with fixed or steady-state initial conditions.
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Fluid/HeatingSystem.png\" border=\"1\"
     alt=\"HeatingSystem.png\">
</p>

<p>
With respect to previous versions, the design
of the connectors has been changed in a non-backward compatible way,
using the recently developed concept
of stream connectors that results in much more reliable simulations
(see also <a href=\"modelica://Modelica/Resources/Documentation/Fluid/Stream-Connectors-Overview-Rationale.pdf\">Stream-Connectors-Overview-Rationale.pdf</a>).
This extension was included in Modelica 3.1.
</p>

<p>
The following parts are useful, when newly starting with this library:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Fluid.UsersGuide\">Modelica.Fluid.UsersGuide</a>.</li>
<li> <a href=\"modelica://Modelica.Fluid.UsersGuide.ReleaseNotes\">Modelica.Fluid.UsersGuide.ReleaseNotes</a>
     summarizes the changes of the library releases.</li>
<li> <a href=\"modelica://Modelica.Fluid.Examples\">Modelica.Fluid.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
Copyright &copy; 2002-2020, Modelica Association and contributors
</p>
</html>"));
  end Fluid;

  package Media "Library of media property models"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;
    import Cv = Modelica.Units.Conversions;

  package Interfaces "Interfaces for media models"
    extends Modelica.Icons.InterfacesPackage;

    partial package PartialMedium
      "Partial medium properties (base package of all media packages)"
      extends Modelica.Media.Interfaces.Types;
      extends Modelica.Icons.MaterialPropertiesPackage;

      // Constants to be set in Medium
      constant Modelica.Media.Interfaces.Choices.IndependentVariables
        ThermoStates "Enumeration type for independent variables";
      constant String mediumName="unusablePartialMedium" "Name of the medium";
      constant String substanceNames[:]={mediumName}
        "Names of the mixture substances. Set substanceNames={mediumName} if only one substance.";
      constant String extraPropertiesNames[:]=fill("", 0)
        "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
      constant Boolean singleState
        "= true, if u and d are not a function of pressure";
      constant Boolean reducedX=true
        "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
      constant Boolean fixedX=false
        "= true if medium contains the equation X = reference_X";
      constant AbsolutePressure reference_p=101325
        "Reference pressure of Medium: default 1 atmosphere";
      constant Temperature reference_T=298.15
        "Reference temperature of Medium: default 25 deg Celsius";
      constant MassFraction reference_X[nX]=fill(1/nX, nX)
        "Default mass fractions of medium";
      constant AbsolutePressure p_default=101325
        "Default value for pressure of medium (for initialization)";
      constant Temperature T_default=Modelica.Units.Conversions.from_degC(20)
        "Default value for temperature of medium (for initialization)";
      constant SpecificEnthalpy h_default=specificEnthalpy_pTX(
              p_default,
              T_default,
              X_default)
        "Default value for specific enthalpy of medium (for initialization)";
      constant MassFraction X_default[nX]=reference_X
        "Default value for mass fractions of medium (for initialization)";
      constant ExtraProperty C_default[nC]=fill(0, nC)
        "Default value for trace substances of medium (for initialization)";

      final constant Integer nS=size(substanceNames, 1) "Number of substances";
      constant Integer nX=nS "Number of mass fractions";
      constant Integer nXi=if fixedX then 0 else if reducedX then nS - 1 else nS
        "Number of structurally independent mass fractions (see docu for details)";

      final constant Integer nC=size(extraPropertiesNames, 1)
        "Number of extra (outside of standard mass-balance) transported properties";
      constant Real C_nominal[nC](min=fill(Modelica.Constants.eps, nC)) = 1.0e-6*
        ones(nC) "Default for the nominal values for the extra properties";
      replaceable record FluidConstants =
          Modelica.Media.Interfaces.Types.Basic.FluidConstants
        "Critical, triple, molecular and other standard data of fluid";

      replaceable record ThermodynamicState
        "Minimal variable set that is available as input argument to every medium function"
        extends Modelica.Icons.Record;
      end ThermodynamicState;

      replaceable partial model BaseProperties
        "Base properties (p, d, T, h, u, R_s, MM and, if applicable, X and Xi) of a medium"
        InputAbsolutePressure p "Absolute pressure of medium";
        InputMassFraction[nXi] Xi(start=reference_X[1:nXi])
          "Structurally independent mass fractions";
        InputSpecificEnthalpy h "Specific enthalpy of medium";
        Density d "Density of medium";
        Temperature T "Temperature of medium";
        MassFraction[nX] X(start=reference_X)
          "Mass fractions (= (component mass)/total mass  m_i/m)";
        SpecificInternalEnergy u "Specific internal energy of medium";
        SpecificHeatCapacity R_s "Gas constant (of mixture if applicable)";
        MolarMass MM "Molar mass (of mixture or single fluid)";
        ThermodynamicState state
          "Thermodynamic state record for optional functions";
        parameter Boolean preferredMediumStates=false
          "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
          annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter Boolean standardOrderComponents=true
          "If true, and reducedX = true, the last element of X will be computed from the other ones";
        Modelica.Units.NonSI.Temperature_degC T_degC=
            Modelica.Units.Conversions.to_degC(T)
          "Temperature of medium in [degC]";
        Modelica.Units.NonSI.Pressure_bar p_bar=
            Modelica.Units.Conversions.to_bar(p)
          "Absolute pressure of medium in [bar]";

        // Local connector definition, used for equation balancing check
        connector InputAbsolutePressure = input SI.AbsolutePressure
          "Pressure as input signal connector";
        connector InputSpecificEnthalpy = input SI.SpecificEnthalpy
          "Specific enthalpy as input signal connector";
        connector InputMassFraction = input SI.MassFraction
          "Mass fraction as input signal connector";

      equation
        if standardOrderComponents then
          Xi = X[1:nXi];

          if fixedX then
            X = reference_X;
          end if;
          if reducedX and not fixedX then
            X[nX] = 1 - sum(Xi);
          end if;
          for i in 1:nX loop
            assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X[" +
              String(i) + "] = " + String(X[i]) + "of substance " +
              substanceNames[i] + "\nof medium " + mediumName +
              " is not in the range 0..1");
          end for;

        end if;

        assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" +
          mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,255}), Text(
                extent={{-152,164},{152,102}},
                textString="%name",
                textColor={0,0,255})}), Documentation(info="<html>
<p>
Model <strong>BaseProperties</strong> is a model within package <strong>PartialMedium</strong>
and contains the <strong>declarations</strong> of the minimum number of
variables that every medium model is supposed to support.
A specific medium inherits from model <strong>BaseProperties</strong> and provides
the equations for the basic properties.</p>
<p>
The BaseProperties model contains the following <strong>7+nXi variables</strong>
(nXi is the number of independent mass fractions defined in package
PartialMedium):
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Variable</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>T</td>
      <td>K</td>
      <td>Temperature</td></tr>
  <tr><td>p</td>
      <td>Pa</td>
      <td>Absolute pressure</td></tr>
  <tr><td>d</td>
      <td>kg/m3</td>
      <td>Density</td></tr>
  <tr><td>h</td>
      <td>J/kg</td>
      <td>Specific enthalpy</td></tr>
  <tr><td>u</td>
      <td>J/kg</td>
      <td>Specific internal energy</td></tr>
  <tr><td>Xi[nXi]</td>
      <td>kg/kg</td>
      <td>Structurally independent mass fractions</td></tr>
  <tr><td>R_s</td>
      <td>J/(kg.K)</td>
      <td>Specific gas constant (of mixture if applicable)</td></tr>
  <tr><td>MM</td>
      <td>kg/mol</td>
      <td>Molar mass</td></tr>
</table>
<p>
In order to implement an actual medium model, one can extend from this
base model and add <strong>5 equations</strong> that provide relations among
these variables. Equations will also have to be added in order to
set all the variables within the ThermodynamicState record state.</p>
<p>
If standardOrderComponents=true, the full composition vector X[nX]
is determined by the equations contained in this base class, depending
on the independent mass fraction vector Xi[nXi].</p>
<p>Additional <strong>2 + nXi</strong> equations will have to be provided
when using the BaseProperties model, in order to fully specify the
thermodynamic conditions. The input connector qualifier applied to
p, h, and nXi indirectly declares the number of missing equations,
permitting advanced equation balance checking by Modelica tools.
Please note that this doesn't mean that the additional equations
should be connection equations, nor that exactly those variables
should be supplied, in order to complete the model.
For further information, see the <a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media User's guide</a>, and
<a href=\"https://specification.modelica.org/v3.4/Ch4.html#balanced-models\">Section 4.7 (Balanced Models) of the Modelica 3.4 specification</a>.</p>
</html>"));
      end BaseProperties;

      replaceable partial function setState_pTX
        "Return thermodynamic state as function of p, T and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_pTX;

      replaceable partial function setState_phX
        "Return thermodynamic state as function of p, h and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_phX;

      replaceable partial function setState_psX
        "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_psX;

      replaceable partial function setState_dTX
        "Return thermodynamic state as function of d, T and composition X or Xi"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_dTX;

      replaceable partial function setSmoothState
        "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;
        input Real x "m_flow or dp";
        input ThermodynamicState state_a "Thermodynamic state if x > 0";
        input ThermodynamicState state_b "Thermodynamic state if x < 0";
        input Real x_small(min=0)
          "Smooth transition in the region -x_small < x < x_small";
        output ThermodynamicState state
          "Smooth thermodynamic state for all x (continuous and differentiable)";
        annotation (Documentation(info="<html>
<p>
This function is used to approximate the equation
</p>
<blockquote><pre>
state = <strong>if</strong> x &gt; 0 <strong>then</strong> state_a <strong>else</strong> state_b;
</pre></blockquote>

<p>
by a smooth characteristic, so that the expression is continuous and differentiable:
</p>

<blockquote><pre>
state := <strong>smooth</strong>(1, <strong>if</strong> x &gt;  x_small <strong>then</strong> state_a <strong>else</strong>
                   <strong>if</strong> x &lt; -x_small <strong>then</strong> state_b <strong>else</strong> f(state_a, state_b));
</pre></blockquote>

<p>
This is performed by applying function <strong>Media.Common.smoothStep</strong>(..)
on every element of the thermodynamic state record.
</p>

<p>
If <strong>mass fractions</strong> X[:] are approximated with this function then this can be performed
for all <strong>nX</strong> mass fractions, instead of applying it for nX-1 mass fractions and computing
the last one by the mass fraction constraint sum(X)=1. The reason is that the approximating function has the
property that sum(state.X) = 1, provided sum(state_a.X) = sum(state_b.X) = 1.
This can be shown by evaluating the approximating function in the abs(x) &lt; x_small
region (otherwise state.X is either state_a.X or state_b.X):
</p>

<blockquote><pre>
X[1]  = smoothStep(x, X_a[1] , X_b[1] , x_small);
X[2]  = smoothStep(x, X_a[2] , X_b[2] , x_small);
   ...
X[nX] = smoothStep(x, X_a[nX], X_b[nX], x_small);
</pre></blockquote>

<p>
or
</p>

<blockquote><pre>
X[1]  = c*(X_a[1]  - X_b[1])  + (X_a[1]  + X_b[1])/2
X[2]  = c*(X_a[2]  - X_b[2])  + (X_a[2]  + X_b[2])/2;
   ...
X[nX] = c*(X_a[nX] - X_b[nX]) + (X_a[nX] + X_b[nX])/2;
c     = (x/x_small)*((x/x_small)^2 - 3)/4
</pre></blockquote>

<p>
Summing all mass fractions together results in
</p>

<blockquote><pre>
sum(X) = c*(sum(X_a) - sum(X_b)) + (sum(X_a) + sum(X_b))/2
       = c*(1 - 1) + (1 + 1)/2
       = 1
</pre></blockquote>

</html>"));
      end setSmoothState;

      replaceable partial function dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DynamicViscosity eta "Dynamic viscosity";
      end dynamicViscosity;

      replaceable partial function thermalConductivity
        "Return thermal conductivity"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output ThermalConductivity lambda "Thermal conductivity";
      end thermalConductivity;

      replaceable function prandtlNumber "Return the Prandtl number"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output PrandtlNumber Pr "Prandtl number";
      algorithm
        Pr := dynamicViscosity(state)*specificHeatCapacityCp(state)/
          thermalConductivity(state);
      end prandtlNumber;

      replaceable partial function pressure "Return pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output AbsolutePressure p "Pressure";
      end pressure;

      replaceable partial function temperature "Return temperature"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Temperature T "Temperature";
      end temperature;

      replaceable partial function density "Return density"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density d "Density";
      end density;

      replaceable partial function specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnthalpy h "Specific enthalpy";
      end specificEnthalpy;

      replaceable partial function specificInternalEnergy
        "Return specific internal energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy u "Specific internal energy";
      end specificInternalEnergy;

      replaceable partial function specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEntropy s "Specific entropy";
      end specificEntropy;

      replaceable partial function specificGibbsEnergy
        "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy g "Specific Gibbs energy";
      end specificGibbsEnergy;

      replaceable partial function specificHelmholtzEnergy
        "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy f "Specific Helmholtz energy";
      end specificHelmholtzEnergy;

      replaceable partial function specificHeatCapacityCp
        "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificHeatCapacity cp
          "Specific heat capacity at constant pressure";
      end specificHeatCapacityCp;

      function heatCapacity_cp = specificHeatCapacityCp
        "Alias for deprecated name";

      replaceable partial function specificHeatCapacityCv
        "Return specific heat capacity at constant volume"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificHeatCapacity cv
          "Specific heat capacity at constant volume";
      end specificHeatCapacityCv;

      function heatCapacity_cv = specificHeatCapacityCv
        "Alias for deprecated name";

      replaceable partial function isentropicExponent
        "Return isentropic exponent"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output IsentropicExponent gamma "Isentropic exponent";
      end isentropicExponent;

      replaceable partial function isentropicEnthalpy
        "Return isentropic enthalpy"
        extends Modelica.Icons.Function;
        input AbsolutePressure p_downstream "Downstream pressure";
        input ThermodynamicState refState "Reference state for entropy";
        output SpecificEnthalpy h_is "Isentropic enthalpy";
        annotation (Documentation(info="<html>
<p>
This function computes an isentropic state transformation:
</p>
<ol>
<li> A medium is in a particular state, refState.</li>
<li> The enthalpy at another state (h_is) shall be computed
     under the assumption that the state transformation from refState to h_is
     is performed with a change of specific entropy ds = 0 and the pressure of state h_is
     is p_downstream and the composition X upstream and downstream is assumed to be the same.</li>
</ol>

</html>"));
      end isentropicEnthalpy;

      replaceable partial function velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output VelocityOfSound a "Velocity of sound";
      end velocityOfSound;

      replaceable partial function isobaricExpansionCoefficient
        "Return overall the isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
        annotation (Documentation(info="<html>
<blockquote><pre>
beta is defined as  1/v * der(v,T), with v = 1/d, at constant pressure p.
</pre></blockquote>
</html>"));
      end isobaricExpansionCoefficient;

      function beta = isobaricExpansionCoefficient
        "Alias for isobaricExpansionCoefficient for user convenience";

      replaceable partial function isothermalCompressibility
        "Return overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SI.IsothermalCompressibility kappa "Isothermal compressibility";
        annotation (Documentation(info="<html>
<blockquote><pre>

kappa is defined as - 1/v * der(v,p), with v = 1/d at constant temperature T.

</pre></blockquote>
</html>"));
      end isothermalCompressibility;

      function kappa = isothermalCompressibility
        "Alias of isothermalCompressibility for user convenience";

      // explicit derivative functions for finite element models
      replaceable partial function density_derp_h
        "Return density derivative w.r.t. pressure at const specific enthalpy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByPressure ddph "Density derivative w.r.t. pressure";
      end density_derp_h;

      replaceable partial function density_derh_p
        "Return density derivative w.r.t. specific enthalpy at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByEnthalpy ddhp
          "Density derivative w.r.t. specific enthalpy";
      end density_derh_p;

      replaceable partial function density_derp_T
        "Return density derivative w.r.t. pressure at const temperature"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByPressure ddpT "Density derivative w.r.t. pressure";
      end density_derp_T;

      replaceable partial function density_derT_p
        "Return density derivative w.r.t. temperature at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByTemperature ddTp
          "Density derivative w.r.t. temperature";
      end density_derT_p;

      replaceable partial function density_derX
        "Return density derivative w.r.t. mass fraction"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density[nX] dddX "Derivative of density w.r.t. mass fraction";
      end density_derX;

      replaceable partial function molarMass
        "Return the molar mass of the medium"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output MolarMass MM "Mixture molar mass";
      end molarMass;

      replaceable function specificEnthalpy_pTX
        "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_pTX(
                p,
                T,
                X));
        annotation (inverse(T=temperature_phX(
                      p,
                      h,
                      X)));
      end specificEnthalpy_pTX;

      replaceable function specificEntropy_pTX
        "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEntropy s "Specific entropy";
      algorithm
        s := specificEntropy(setState_pTX(
                p,
                T,
                X));

        annotation (inverse(T=temperature_psX(
                      p,
                      s,
                      X)));
      end specificEntropy_pTX;

      replaceable function density_pTX "Return density from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:] "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                X));
      end density_pTX;

      replaceable function temperature_phX
        "Return temperature from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_phX(
                p,
                h,
                X));
      end temperature_phX;

      replaceable function density_phX "Return density from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_phX(
                p,
                h,
                X));
      end density_phX;

      replaceable function temperature_psX
        "Return temperature from p,s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_psX(
                p,
                s,
                X));
        annotation (inverse(s=specificEntropy_pTX(
                      p,
                      T,
                      X)));
      end temperature_psX;

      replaceable function density_psX "Return density from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_psX(
                p,
                s,
                X));
      end density_psX;

      replaceable function specificEnthalpy_psX
        "Return specific enthalpy from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_psX(
                p,
                s,
                X));
      end specificEnthalpy_psX;

      type MassFlowRate = SI.MassFlowRate (
          quantity="MassFlowRate." + mediumName,
          min=-1.0e5,
          max=1.e5) "Type for mass flow rate with medium specific attributes";

      annotation (Documentation(info="<html>
<p>
<strong>PartialMedium</strong> is a package and contains all <strong>declarations</strong> for
a medium. This means that constants, models, and functions
are defined that every medium is supposed to support
(some of them are optional). A medium package
inherits from <strong>PartialMedium</strong> and provides the
equations for the medium. The details of this package
are described in
<a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.
</p>
</html>",   revisions="<html>

</html>"));
    end PartialMedium;

    package Choices "Types, constants to define menu choices"
      extends Modelica.Icons.Package;

      type IndependentVariables = enumeration(
          T "Temperature",
          pT "Pressure, Temperature",
          ph "Pressure, Specific Enthalpy",
          phX "Pressure, Specific Enthalpy, Mass Fraction",
          pTX "Pressure, Temperature, Mass Fractions",
          dTX "Density, Temperature, Mass Fractions")
        "Enumeration defining the independent variables of a medium";
      annotation (Documentation(info="<html>
<p>
Enumerations and data types for all types of fluids
</p>

<p>
Note: Reference enthalpy might have to be extended with enthalpy of formation.
</p>
</html>"));
    end Choices;

    package Types "Types to be used in fluid models"
      extends Modelica.Icons.Package;

      type AbsolutePressure = SI.AbsolutePressure (
          min=0,
          max=1.e8,
          nominal=1.e5,
          start=1.e5)
        "Type for absolute pressure with medium specific attributes";

      type Density = SI.Density (
          min=0,
          max=1.e5,
          nominal=1,
          start=1) "Type for density with medium specific attributes";

      type DynamicViscosity = SI.DynamicViscosity (
          min=0,
          max=1.e8,
          nominal=1.e-3,
          start=1.e-3)
        "Type for dynamic viscosity with medium specific attributes";

      type MassFraction = Real (
          quantity="MassFraction",
          final unit="kg/kg",
          min=0,
          max=1,
          nominal=0.1) "Type for mass fraction with medium specific attributes";

      type MolarMass = SI.MolarMass (
          min=0.001,
          max=0.25,
          nominal=0.032) "Type for molar mass with medium specific attributes";

      type IsentropicExponent = SI.RatioOfSpecificHeatCapacities (
          min=1,
          max=500000,
          nominal=1.2,
          start=1.2)
        "Type for isentropic exponent with medium specific attributes";

      type SpecificEnergy = SI.SpecificEnergy (
          min=-1.0e8,
          max=1.e8,
          nominal=1.e6)
        "Type for specific energy with medium specific attributes";

      type SpecificInternalEnergy = SpecificEnergy
        "Type for specific internal energy with medium specific attributes";

      type SpecificEnthalpy = SI.SpecificEnthalpy (
          min=-1.0e10,
          max=1.e10,
          nominal=1.e6)
        "Type for specific enthalpy with medium specific attributes";

      type SpecificEntropy = SI.SpecificEntropy (
          min=-1.e7,
          max=1.e7,
          nominal=1.e3)
        "Type for specific entropy with medium specific attributes";

      type SpecificHeatCapacity = SI.SpecificHeatCapacity (
          min=0,
          max=1.e7,
          nominal=1.e3,
          start=1.e3)
        "Type for specific heat capacity with medium specific attributes";

      type Temperature = SI.Temperature (
          min=1,
          max=1.e4,
          nominal=300,
          start=288.15) "Type for temperature with medium specific attributes";

      type ThermalConductivity = SI.ThermalConductivity (
          min=0,
          max=500,
          nominal=1,
          start=1)
        "Type for thermal conductivity with medium specific attributes";

      type PrandtlNumber = SI.PrandtlNumber (
          min=1e-3,
          max=1e5,
          nominal=1.0) "Type for Prandtl number with medium specific attributes";

      type VelocityOfSound = SI.Velocity (
          min=0,
          max=1.e5,
          nominal=1000,
          start=1000)
        "Type for velocity of sound with medium specific attributes";

      type ExtraProperty = Real (min=0.0, start=1.0)
        "Type for unspecified, mass-specific property transported by flow";

      type IsobaricExpansionCoefficient = Real (
          min=0,
          max=1.0e8,
          unit="1/K")
        "Type for isobaric expansion coefficient with medium specific attributes";

      type DerDensityByPressure = SI.DerDensityByPressure
        "Type for partial derivative of density with respect to pressure with medium specific attributes";

      type DerDensityByEnthalpy = SI.DerDensityByEnthalpy
        "Type for partial derivative of density with respect to enthalpy with medium specific attributes";

      type DerDensityByTemperature = SI.DerDensityByTemperature
        "Type for partial derivative of density with respect to temperature with medium specific attributes";

      package Basic
      "The most basic version of a record used in several degrees of detail"
        extends Icons.Package;

        record FluidConstants
          "Critical, triple, molecular and other standard data of fluid"
          extends Modelica.Icons.Record;
          String iupacName
            "Complete IUPAC name (or common name, if non-existent)";
          String casRegistryNumber
            "Chemical abstracts sequencing number (if it exists)";
          String chemicalFormula
            "Chemical formula, (brutto, nomenclature according to Hill";
          String structureFormula "Chemical structure formula";
          MolarMass molarMass "Molar mass";
        end FluidConstants;
      end Basic;
    end Types;
    annotation (Documentation(info="<html>
<p>
This package provides basic interfaces definitions of media models for different
kind of media.
</p>
</html>"));
  end Interfaces;
  annotation (preferredView="info",Documentation(info="<html>
<p>
This library contains <a href=\"modelica://Modelica.Media.Interfaces\">interface</a>
definitions for media and the following <strong>property</strong> models for
single and multiple substance fluids with one and multiple phases:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Media.IdealGases\">Ideal gases:</a><br>
     1241 high precision gas models based on the
     NASA Glenn coefficients, plus ideal gas mixture models based
     on the same data.</li>
<li> <a href=\"modelica://Modelica.Media.Water\">Water models:</a><br>
     ConstantPropertyLiquidWater, WaterIF97 (high precision
     water model according to the IAPWS/IF97 standard)</li>
<li> <a href=\"modelica://Modelica.Media.Air\">Air models:</a><br>
     SimpleAir, DryAirNasa, ReferenceAir, MoistAir, ReferenceMoistAir.</li>
<li> <a href=\"modelica://Modelica.Media.Incompressible\">
     Incompressible media:</a><br>
     TableBased incompressible fluid models (properties are defined by tables rho(T),
     HeatCapacity_cp(T), etc.)</li>
<li> <a href=\"modelica://Modelica.Media.CompressibleLiquids\">
     Compressible liquids:</a><br>
     Simple liquid models with linear compressibility</li>
<li> <a href=\"modelica://Modelica.Media.R134a\">Refrigerant Tetrafluoroethane (R134a)</a>.</li>
</ul>
<p>
The following parts are useful, when newly starting with this library:</p>
<ul>
<li> <a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.MediumUsage\">Modelica.Media.UsersGuide.MediumUsage</a>
     describes how to use a medium model in a component model.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.MediumDefinition\">
     Modelica.Media.UsersGuide.MediumDefinition</a>
     describes how a new fluid medium model has to be implemented.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.ReleaseNotes\">Modelica.Media.UsersGuide.ReleaseNotes</a>
     summarizes the changes of the library releases.</li>
<li> <a href=\"modelica://Modelica.Media.Examples\">Modelica.Media.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>February 01, 2017</em> by Thomas Beutlich:<br>
    Fixed data errors of the NASA Glenn coefficients in some ideal gases (CH2, CH3, CH3OOH, C2CL2, C2CL4, C2CL6, C2HCL, C2HCL3, CH2CO_ketene, O_CH_2O, HO_CO_2OH, CH2BrminusCOOH, C2H3CL, CH2CLminusCOOH, HO2, HO2minus, OD, ODminus), see <a href=\"https://github.com/modelica/ModelicaStandardLibrary/issues/1922\">#1922</a></li>
<li><em>May 16, 2013</em> by Stefan Wischhusen (XRG Simulation):<br>
    Added new media models Air.ReferenceMoistAir, Air.ReferenceAir, R134a.</li>
<li><em>May 25, 2011</em> by Francesco Casella:<br>Added min/max attributes to Water, TableBased, MixtureGasNasa, SimpleAir and MoistAir local types.</li>
<li><em>May 25, 2011</em> by Stefan Wischhusen:<br>Added individual settings for polynomial fittings of properties.</li>
</ul>
</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics={
          Line(
            points = {{-76,-80},{-62,-30},{-32,40},{4,66},{48,66},{73,45},{62,-8},{48,-50},{38,-80}},
            color={64,64,64},
            smooth=Smooth.Bezier),
          Line(
            points={{-40,20},{68,20}},
            color={175,175,175}),
          Line(
            points={{-40,20},{-44,88},{-44,88}},
            color={175,175,175}),
          Line(
            points={{68,20},{86,-58}},
            color={175,175,175}),
          Line(
            points={{-60,-28},{56,-28}},
            color={175,175,175}),
          Line(
            points={{-60,-28},{-74,84},{-74,84}},
            color={175,175,175}),
          Line(
            points={{56,-28},{70,-80}},
            color={175,175,175}),
          Line(
            points={{-76,-80},{38,-80}},
            color={175,175,175}),
          Line(
            points={{-76,-80},{-94,-16},{-94,-16}},
            color={175,175,175})}));
  end Media;

  package Thermal
  "Library of thermal system components to model heat transfer and simple thermo-fluid pipe flow"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;

    package HeatTransfer
    "Library of 1-dimensional heat transfer with lumped elements"
      extends Modelica.Icons.Package;

      package Interfaces "Connectors and partial models"
        extends Modelica.Icons.InterfacesPackage;

        partial connector HeatPort "Thermal port for 1-dim. heat transfer"
          SI.Temperature T "Port temperature";
          flow SI.HeatFlowRate Q_flow
            "Heat flow rate (positive if flowing from outside into the component)";
          annotation (Documentation(info="<html>

</html>"));
        end HeatPort;

        connector HeatPort_a
          "Thermal port for 1-dim. heat transfer (filled rectangular icon)"

          extends HeatPort;

          annotation(defaultComponentName = "port_a",
            Documentation(info="<html>
<p>This connector is used for 1-dimensional heat flow between components.
The variables in the connector are:</p>
<blockquote><pre>
T       Temperature in [Kelvin].
Q_flow  Heat flow rate in [Watt].
</pre></blockquote>
<p>According to the Modelica sign convention, a <strong>positive</strong> heat flow
rate <strong>Q_flow</strong> is considered to flow <strong>into</strong> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <strong>HeatPort_a</strong> and
<strong>HeatPort_b</strong> are identical with the only exception of the different
<strong>icon layout</strong>.</p></html>"),     Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics={Rectangle(
                  extent={{-50,50},{50,-50}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-120,120},{100,60}},
                  textColor={191,0,0},
                  textString="%name")}));
        end HeatPort_a;
        annotation (Documentation(info="<html>

</html>"));
      end Interfaces;
      annotation (
         Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}), graphics={
          Polygon(
            origin = {13.758,27.517},
            lineColor = {128,128,128},
            fillColor = {192,192,192},
            fillPattern = FillPattern.Solid,
            points = {{-54,-6},{-61,-7},{-75,-15},{-79,-24},{-80,-34},{-78,-42},{-73,-49},{-64,-51},{-57,-51},{-47,-50},{-41,-43},{-38,-35},{-40,-27},{-40,-20},{-42,-13},{-47,-7},{-54,-5},{-54,-6}}),
        Polygon(
            origin = {13.758,27.517},
            fillColor = {160,160,164},
            fillPattern = FillPattern.Solid,
            points = {{-75,-15},{-79,-25},{-80,-34},{-78,-42},{-72,-49},{-64,-51},{-57,-51},{-47,-50},{-57,-47},{-65,-45},{-71,-40},{-74,-33},{-76,-23},{-75,-15},{-75,-15}}),
          Polygon(
            origin = {13.758,27.517},
            lineColor = {160,160,164},
            fillColor = {192,192,192},
            fillPattern = FillPattern.Solid,
            points = {{39,-6},{32,-7},{18,-15},{14,-24},{13,-34},{15,-42},{20,-49},{29,-51},{36,-51},{46,-50},{52,-43},{55,-35},{53,-27},{53,-20},{51,-13},{46,-7},{39,-5},{39,-6}}),
          Polygon(
            origin = {13.758,27.517},
            fillColor = {160,160,164},
            fillPattern = FillPattern.Solid,
            points = {{18,-15},{14,-25},{13,-34},{15,-42},{21,-49},{29,-51},{36,-51},{46,-50},{36,-47},{28,-45},{22,-40},{19,-33},{17,-23},{18,-15},{18,-15}}),
          Polygon(
            origin = {13.758,27.517},
            lineColor = {191,0,0},
            fillColor = {191,0,0},
            fillPattern = FillPattern.Solid,
            points = {{-9,-23},{-9,-10},{18,-17},{-9,-23}}),
          Line(
            origin = {13.758,27.517},
            points = {{-41,-17},{-9,-17}},
            color = {191,0,0},
            thickness = 0.5),
          Line(
            origin = {13.758,27.517},
            points = {{-17,-40},{15,-40}},
            color = {191,0,0},
            thickness = 0.5),
          Polygon(
            origin = {13.758,27.517},
            lineColor = {191,0,0},
            fillColor = {191,0,0},
            fillPattern = FillPattern.Solid,
            points = {{-17,-46},{-17,-34},{-40,-40},{-17,-46}})}),
                                Documentation(info="<html>
<p>
This package contains components to model <strong>1-dimensional heat transfer</strong>
with lumped elements.</p>
</html>",     revisions="<html>
<ul>
<li><em>July 15, 2002</em>
       by Michael Tiller, <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Nikolaus Sch&uuml;rmann:<br>
       Implemented.
</li>
<li><em>June 13, 2005</em>
       by <a href=\"https://www.haumer.at/\">Anton Haumer</a><br>
       Refined placing of connectors (cosmetic).<br>
       Refined all Examples; removed Examples.FrequencyInverter, introducing Examples.Motor<br>
       Introduced temperature dependent correction (1 + alpha*(T - T_ref)) in Fixed/PrescribedHeatFlow<br>
</li>
  <li> v1.1.1 2007/11/13 Anton Haumer<br>
       components moved to sub-packages</li>
  <li> v1.2.0 2009/08/26 Anton Haumer<br>
       added component ThermalCollector</li>

</ul>
</html>"));
    end HeatTransfer;
    annotation (
     Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
      Line(
      origin={-47.5,11.6667},
      points={{-2.5,-91.6667},{17.5,-71.6667},{-22.5,-51.6667},{17.5,-31.6667},{-22.5,-11.667},{17.5,8.3333},{-2.5,28.3333},{-2.5,48.3333}},
        smooth=Smooth.Bezier),
      Polygon(
      origin={-50.0,68.333},
      pattern=LinePattern.None,
      fillPattern=FillPattern.Solid,
        points={{0.0,21.667},{-10.0,-8.333},{10.0,-8.333}}),
      Line(
      origin={2.5,11.6667},
      points={{-2.5,-91.6667},{17.5,-71.6667},{-22.5,-51.6667},{17.5,-31.6667},{-22.5,-11.667},{17.5,8.3333},{-2.5,28.3333},{-2.5,48.3333}},
        smooth=Smooth.Bezier),
      Polygon(
      origin={0.0,68.333},
      pattern=LinePattern.None,
      fillPattern=FillPattern.Solid,
        points={{0.0,21.667},{-10.0,-8.333},{10.0,-8.333}}),
      Line(
      origin={52.5,11.6667},
      points={{-2.5,-91.6667},{17.5,-71.6667},{-22.5,-51.6667},{17.5,-31.6667},{-22.5,-11.667},{17.5,8.3333},{-2.5,28.3333},{-2.5,48.3333}},
        smooth=Smooth.Bezier),
      Polygon(
      origin={50.0,68.333},
      pattern=LinePattern.None,
      fillPattern=FillPattern.Solid,
        points={{0.0,21.667},{-10.0,-8.333},{10.0,-8.333}})}),
      Documentation(info="<html>
<p>
This package contains libraries to model heat transfer
and fluid heat flow.
</p>
</html>"));
  end Thermal;

  package Math
  "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices"
    extends Modelica.Icons.Package;

    package Nonlinear "Library of functions operating on nonlinear equations"
      extends Modelica.Icons.Package;

      package Interfaces "Interfaces for functions"
        extends Modelica.Icons.InterfacesPackage;

      encapsulated partial function partialScalarFunction
          "Interface for a function with one input and one output Real signal"
        import Modelica;
        extends Modelica.Icons.Function;
        input Real u "Independent variable";
        output Real y "Dependent variable y=f(u)";
          annotation (Documentation(info="<html>
<p>
This partial function defines the interface of a function with
one input and one output Real signal. The scalar functions
of <a href=\"modelica://Modelica.Math.Nonlinear\">Modelica.Math.Nonlinear</a>
are derived from this base type by inheritance.
This allows to use these functions directly as function arguments
to a function, see, .e.g.,
<a href=\"modelica://Modelica.Math.Nonlinear.Examples\">Math.Nonlinear.Examples</a>.
</p>

</html>"));
      end partialScalarFunction;
        annotation (Documentation(info="<html>
<p>
Interface definitions of functions. The main purpose is to use functions
derived from these interface definitions as function arguments
to a function, see, .e.g.,
<a href=\"modelica://Modelica.Math.Nonlinear.Examples\">Math.Nonlinear.Examples</a>.
</p>
</html>"));
      end Interfaces;
      annotation (Documentation(info="<html>
<p>
This package contains functions to perform tasks such as numerically integrating
a function, or solving a nonlinear algebraic equation system.
The common feature of the functions in this package is
that the nonlinear characteristics are passed as user definable
functions.
</p>

<p>
For details about how to define and to use functions as input arguments
to functions, see
<a href=\"modelica://ModelicaReference.Classes.'function'\">ModelicaReference.Classes.'function'</a>
or <a href=\"https://specification.modelica.org/v3.4/Ch12.html#functional-input-arguments-to-functions\">Section 12.4.2
(Functional Input Arguments to Functions) of the Modelica 3.4 specification</a>.
</p>

</html>",     revisions="<html>
<ul>
<li><em>July 2010 </em> by Martin Otter (DLR-RM):<br>
    Included in MSL 3.2, adapted, and documentation improved</li>

<li><em>March 2010 </em> by Andreas Pfeiffer (DLR-RM):<br>
    Adapted the quadrature function from Gerhard Schillhuber and
    the solution of one non-linear equation in one unknown from
    Modelica.Media.Common.OneNonLinearEquation so that
    function objects are used.</li>

<li><em>June 2002 </em> by Gerhard Schillhuber (master thesis at DLR-RM):<br>
       Adaptive quadrature to compute the curve length of a Spline.</li>
</ul>
</html>"),     Icon(graphics={Polygon(points={{-44,-52},{-44,-26},{-17.1,
                  44.4},{-11.4,52.6},{-5.8,57.1},{-0.2,57.8},{5.4,54.6},{11.1,47.7},
                  {16.7,37.4},{23.1,22.1},{31.17,-0.8},{48,-52},{-44,-52}},
              lineColor={135,135,135},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid)}));
    end Nonlinear;

    package Random "Library of functions for generating random numbers"
       extends Modelica.Icons.Package;

      package Generators
      "Library of functions generating uniform random numbers in the range 0 < random <= 1.0 (with exposed state vectors)"
        extends Modelica.Icons.Package;

        package Xorshift64star "Random number generator xorshift64*"

          constant Integer nState=2 "The dimension of the internal state vector";
          extends Modelica.Icons.Package;

          function initialState
            "Returns an initial state for the xorshift64* algorithm"
            extends Modelica.Icons.Function;
            input Integer localSeed
              "The local seed to be used for generating initial states";
            input Integer globalSeed
              "The global seed to be combined with the local seed";
            output Integer state[nState] "The generated initial states";
          protected
            Real r "Random number not used outside the function";

            /* According to http://vigna.di.unimi.it/ftp/papers/xorshift.pdf, the xorshoft*
     random number generator generates statistically random numbers from a bad seed
      within one iteration. To be on the safe side, 10 iterations are actually used
    */
            constant Integer p = 10 "The number of iterations to use";

          algorithm
            // If seed=0 use a large prime number as seed (seed must be different from 0).
            if localSeed == 0 and globalSeed == 0 then
              state := {126247697,globalSeed};
            else
              state := {localSeed,globalSeed};
            end if;

            // Generate p-times a random number, in order to get a "good" state
            // even if starting from a bad seed.
            for i in 1:p loop
              (r,state) := random(state);
            end for;
          annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
state = Xorshift64star.<strong>initialState</strong>(localSeed, globalSeed);
</pre></blockquote>

<h4>Description</h4>
<p>
Generates the initial state vector <strong>state</strong> for the Xorshift64star random number generator
(= xorshift64* algorithm), from
two Integer numbers given as input (arguments localSeed, globalSeed). Any Integer numbers
can be given (including zero or negative number). The function returns
a reasonable initial state vector with the following strategy:
</p>

<p>
If both input
arguments are zero, a fixed non-zero value is used internally for localSeed.
According to <a href=\"http://vigna.di.unimi.it/ftp/papers/xorshift.pdf\">xorshift.pdf</a>,
the xorshift64* random number generator generates statistically random numbers from a
bad seed within one iteration. To be on the safe side, actually 10 random numbers are generated
and the returned state is the one from the last iteration.
</p>

<h4>Example</h4>
<blockquote><pre>
  <strong>parameter</strong> Integer localSeed;
  <strong>parameter</strong> Integer globalSeed;
  Integer state[Xorshift64star.nState];
<strong>initial equation</strong>
  state = initialState(localSeed, globalSeed);
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Generators.Xorshift64star.random\">Random.Generators.Xorshift64star.random</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
          end initialState;

          pure function random
            "Returns a uniform random number with the xorshift64* algorithm"
            extends Modelica.Icons.Function;
            input Integer stateIn[nState]
              "The internal states for the random number generator";
            output Real result
              "A random number with a uniform distribution on the interval (0,1]";
            output Integer stateOut[nState]
              "The new internal states of the random number generator";
            external "C" ModelicaRandom_xorshift64star(stateIn, stateOut, result)
              annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaRandom.h\"", Library="ModelicaExternalC");
            annotation(Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(r, stateOut) = Xorshift64star.<strong>random</strong>(stateIn);
</pre></blockquote>

<h4>Description</h4>
<p>
Returns a uniform random number r in the range 0 &lt; r &le; 1 with the xorshift64* algorithm.
Input argument <strong>stateIn</strong> is the state vector of the previous call.
Output argument <strong>stateOut</strong> is the updated state vector.
If the function is called with identical stateIn vectors, exactly the
same random number r is returned.
</p>

<h4>Example</h4>
<blockquote><pre>
  <strong>parameter</strong> Integer localSeed;
  <strong>parameter</strong> Integer globalSeed;
  Real r;
  Integer state[Xorshift64star.nState];
<strong>initial equation</strong>
  state = initialState(localSeed, globalSeed);
<strong>equation</strong>
  <strong>when</strong> sample(0,0.1) <strong>then</strong>
    (r, state) = random(<strong>pre</strong>(state));
  <strong>end when</strong>;
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Generators.Xorshift64star.initialState\">Random.Generators.Xorshift64star.initialState</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
          end random;
          annotation (Documentation(info="<html>
<p>
Random number generator <strong>xorshift64*</strong>. This generator has a period of 2^64
(the period defines the number of random numbers generated before the sequence begins to repeat itself).
For an overview, comparison with other random number generators, and links to articles, see
<a href=\"modelica://Modelica.Math.Random.Generators\">Math.Random.Generators</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"),     Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={
            Ellipse(
              extent={{-64,0},{-14,-50}},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{12,52},{62,2}},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid)}));
        end Xorshift64star;

        package Xorshift128plus "Random number generator xorshift128+"

          constant Integer nState=4 "The dimension of the internal state vector";
          extends Modelica.Icons.Package;

          function initialState
            "Returns an initial state for the xorshift128+ algorithm"
            extends Modelica.Icons.Function;
            input Integer localSeed
              "The local seed to be used for generating initial states";
            input Integer globalSeed
              "The global seed to be combined with the local seed";
            output Integer state[nState] "The generated initial states";
          algorithm
            state := Utilities.initialStateWithXorshift64star(
                    localSeed,
                    globalSeed,
                    size(state, 1));
            annotation(Inline=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
state = Xorshift128plus.<strong>initialState</strong>(localSeed, globalSeed);
</pre></blockquote>

<h4>Description</h4>
<p>
Generates an initial state vector for the Xorshift128plus random number generator
(= xorshift128+ algorithm), from
two Integer numbers given as input (arguments localSeed, globalSeed). Any Integer numbers
can be given (including zero or negative number). The function returns
a reasonable initial state vector with the following strategy:
</p>

<p>
The <a href=\"modelica://Modelica.Math.Random.Generators.Xorshift64star\">Xorshift64star</a>
random number generator is used to fill the internal state vector with 64 bit random numbers.
</p>

<h4>Example</h4>
<blockquote><pre>
  <strong>parameter</strong> Integer localSeed;
  <strong>parameter</strong> Integer globalSeed;
  Integer state[Xorshift128plus.nState];
<strong>initial equation</strong>
  state = initialState(localSeed, globalSeed);
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Generators.Xorshift128plus.random\">Random.Generators.Xorshift128plus.random</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
          end initialState;

          pure function random
            "Returns a uniform random number with the xorshift128+ algorithm"
            extends Modelica.Icons.Function;
            input Integer stateIn[nState]
              "The internal states for the random number generator";
            output Real result
              "A random number with a uniform distribution on the interval (0,1]";
            output Integer stateOut[nState]
              "The new internal states of the random number generator";
            external "C" ModelicaRandom_xorshift128plus(stateIn, stateOut, result)
              annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaRandom.h\"", Library="ModelicaExternalC");
            annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(r, stateOut) = Xorshift128plus.<strong>random</strong>(stateIn);
</pre></blockquote>

<h4>Description</h4>
<p>
Returns a uniform random number in the range 0 &lt; random &le; 1 with the xorshift128+ algorithm.
Input argument <strong>stateIn</strong> is the state vector of the previous call.
Output argument <strong>stateOut</strong> is the updated state vector.
If the function is called with identical stateIn vectors, exactly the
same random number r is returned.
</p>

<h4>Example</h4>
<blockquote><pre>
  <strong>parameter</strong> Integer localSeed;
  <strong>parameter</strong> Integer globalSeed;
  Real r;
  Integer state[Xorshift128plus.nState];
<strong>initial equation</strong>
  state = initialState(localSeed, globalSeed);
<strong>equation</strong>
  <strong>when</strong> sample(0,0.1) <strong>then</strong>
    (r, state) = random(<strong>pre</strong>(state));
  <strong>end when</strong>;
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Generators.Xorshift128plus.initialState\">Random.Generators.Xorshift128plus.initialState</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
          end random;
          annotation (Documentation(info="<html>
<p>
Random number generator <strong>xorshift128+</strong>. This generator has a period of 2^128
(the period defines the number of random numbers generated before the sequence begins to repeat itself).
For an overview, comparison with
other random number generators, and links to articles, see
<a href=\"modelica://Modelica.Math.Random.Generators\">Math.Random.Generators</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"),
         Icon(graphics={
            Ellipse(
              extent={{-70,60},{-20,10}},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{32,58},{82,8}},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-20,-12},{30,-62}},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid)}));
        end Xorshift128plus;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
              -100,-100},{100,100}}), graphics={Line(
            points={{-90,-54},{-50,-54},{-50,54},{50,54},{50,-54},{84,-54}})}), Documentation(info="<html>
<p>
This package contains various pseudo random number generators. A random number generator is a package
that consists of the following elements:
</p>
<ul>
<li> Integer <strong>nState</strong> is a constant that defines the length of the internal state vector
     (in order that an appropriate Integer vector of this length can be declared, depending on
     the selected random number generator).</li>
<li> Function <strong>initialState(..)</strong> is used to initialize the state of the random number generator
     by providing Integer seeds and calling the random number generator often enough that
     statistically relevant random numbers are returned by every call of function random(..).</li>
<li> Function <strong>random(..)</strong> is used to return a random number of type Real in the range
     0.0 &lt; random &le; 1.0 for every call.
     Furthermore, the updated (internal) state of the random number generator is returned as well.
    </li>
</ul>

<p>
The Generators package contains the <strong>xorshift</strong> suite of random number generators
from Sebastiano Vigna (from 2014; based on work of George Marsaglia).
The properties of these random
number generators are summarized below and compared with the often used
Mersenne Twister (MT19937-64) generator. The table is based on
<a href=\"http://xorshift.di.unimi.it/\">http://xorshift.di.unimi.it/</a> and on the
articles:
</p>
<blockquote>
<p>
Sebastiano Vigna:
<a href=\"http://vigna.di.unimi.it/ftp/papers/xorshift.pdf\">An experimental exploration of Marsaglia's xorshift generators, scrambled</a>, 2014.<br>
Sebastiano Vigna:
<a href=\"http://vigna.di.unimi.it/ftp/papers/xorshiftplus.pdf\">Further scramblings of Marsaglia's xorshift generators</a>, 2014.<br>
</p>
</blockquote>

<p>
Summary of the properties of the random number generators:
</p>

<blockquote>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Property</th>
    <th>xorshift64*</th>
    <th>xorshift128+</th>
    <th>xorshift1024*</th>
    <th>MT19937-64</th></tr>

<tr><td>Period</td>
    <td>2^64</td>
    <td>2^128</td>
    <td>2^1024</td>
    <td>2^19937</td></tr>

<tr><td>Length of state (# 32 bit integer)</td>
    <td>2</td>
    <td>4</td>
    <td>33</td>
    <td>624</td></tr>

<tr><td>Statistic failures (Big Crush)</td>
    <td>363</td>
    <td>64</td>
    <td>51</td>
    <td>516</td></tr>

<tr><td>Systematic failures (Big Crush)</td>
    <td>yes</td>
    <td>no</td>
    <td>no</td>
    <td>yes</td></tr>

<tr><td>Worst case startup</td>
    <td> &gt; 1 call       </td>
    <td> &gt; 20 calls     </td>
    <td> &gt; 100 calls    </td>
    <td> &gt; 100000 calls </td></tr>

<tr><td>Run time (MT=1.0)</td>
    <td> 0.39 </td>
    <td> 0.27 </td>
    <td> 0.33 </td>
    <td> 1.0  </td></tr>
</table>
</blockquote>

<p>
Further explanations of the properties above:
</p>

<ul>
<li> The <strong>period</strong> defines the number of random numbers generated
     before the sequence begins to repeat itself. According to
     \"<a href=\"http://xorshift.di.unimi.it/\">A long period does not imply high quality</a>\"
     a period of 2^1024 is by far large enough for even massively parallel simulations
     with huge number of random number computations per simulation.
     A period of 2^128 might be not enough for massively parallel simulations.
     </li>

<li> <strong>Length of state (# 32 bit integer)</strong> defines the number of \"int\" (that is Modelica Integer) elements
     used for the internal state vector.</li>

<li> <strong>Big Crush</strong> is part of <a href=\"http://simul.iro.umontreal.ca/testu01/tu01.html\">TestU01</a>
     a huge framework for testing random number generators.
     According to these tests, the statistical properties of the xorshift random number
     generators are better than the ones of the Mersenne Twister random number generator.</li>

<li> <strong>Worst case startup</strong> means how many calls are needed until getting
     from a bad seed to random numbers with appropriate statistical properties.
     Here, the xorshift random number suite has much better properties
     than the Mersenne Twister. When initializing a random number generator, the above property
     is taken into account and appropriate random numbers are generated, so that a subsequent
     call of random(..) will generate statistically relevant random numbers, even if the user
     provides a bad initial seed (such as localSeed=1). This means, any Integer number can be given as
     initial seed without influencing the quality of the generated random numbers.</li>

<li> <strong>Run time</strong> shows that the xorshift random number generators are
     all much faster than the Mersenne Twister random number generator, although
     this is not really relevant for most simulations, because the execution
     time of the other parts of the simulations is usually much larger.</li>
</ul>

<p>
The xorshift random number generators are used in the following way in the
<a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a> package:
</p>
<ol>
<li> Xorshift64star (xorshift64*) is used to generate the initial internal state vectors of the
     other generators from two Integer values, due
     to the very good startup properties.</li>

<li> Xorshift128plus (xorshift128+) is the random number generator
     used by the blocks in <a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a>.
     Since these blocks hold the internal state vector for every block instance, and the
     internal state vector is copied whenever a new random number is drawn, it is important
     that the internal state vector is short (and still has good statistical properties
     as shown in the table above).</li>

<li> Xorshift1024star (xorshift1024*) is the basis of the impure function
     <a href=\"modelica://Modelica.Math.Random.Utilities.impureRandom\">Math.Random.Utilities.impureRandom</a>
     which in turn is used with
     <a href=\"modelica://Modelica.Blocks.Noise.GlobalSeed\">Blocks.Noise.GlobalSeed</a>.
     The internal state vector is not exposed. It is updated internally, whenever a new random number
     is drawn.</li>
</ol>

<p>
Note, the generators produce 64 bit random numbers.
These numbers are mapped to the 52 bit mantissa of double numbers in the range 0.0 .. 1.0.
</p>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end Generators;

      package Utilities
      "Library of utility functions for the Random package (usually of no interest for the user)"
        extends Modelica.Icons.UtilitiesPackage;

        function initialStateWithXorshift64star
          "Return an initial state vector for a random number generator (based on xorshift64star algorithm)"
          import Modelica.Math.Random.Generators.Xorshift64star;
          extends Modelica.Icons.Function;
          input Integer localSeed
            "The local seed to be used for generating initial states";
          input Integer globalSeed
            "The global seed to be combined with the local seed";
          input Integer nState(min=1) "The dimension of the state vector (>= 1)";
          output Integer[nState] state "The generated initial states";

        protected
          Real r "Random number only used inside function";
          Integer aux[2] "Intermediate container of state integers";
          Integer nStateEven "Highest even number <= nState";
        algorithm
          // Set the first 2 states by using the initialState() function
          aux            := Xorshift64star.initialState(localSeed, globalSeed);
          if nState >= 2 then
            state[1:2]   := aux;
          else
            state[1]     := aux[1];
          end if;

          // Fill the next elements of the state vector
          nStateEven     := 2*div(nState, 2);
          for i in 3:2:nStateEven loop
            (r,aux)      := Xorshift64star.random(state[i-2:i-1]);
            state[i:i+1] := aux;
          end for;

          // If nState is uneven, fill the last element as well
          if nState >= 3 and nState <> nStateEven then
            (r,aux)       := Xorshift64star.random(state[nState-2:nState-1]);
            state[nState] := aux[1];
          end if;

          annotation (Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>",         info="<html>
<h4>Syntax</h4>
<blockquote><pre>
state = Utilities.<strong>initialStateWithXorshift6star</strong>(localSeed, globalSeed, nState);
</pre></blockquote>

<h4>Description</h4>

<p>
The <a href=\"modelica://Modelica.Math.Random.Generators.Xorshift64star\">Xorshift64star</a>
random number generator is used to fill a state vector of length nState (nState &ge; 1) with random numbers and return
this vector. Arguments localSeed and globalSeed are any Integer numbers (including zero or negative number)
that characterize the initial state.
If the same localSeed, globalSeed, nState is given, the same state vector is returned.
</p>

<h4>Example</h4>
<blockquote><pre>
  parameter Integer localSeed;
  parameter Integer globalSeed;
  Integer state[33];
<strong>initial equation</strong>
  state = Utilities.initialStateWithXorshift64star(localSeed, globalSeed, size(state,1));
</pre></blockquote>
</html>"));
        end initialStateWithXorshift64star;

        impure function automaticGlobalSeed
          "Creates an automatic integer seed (typically from the current time and process id; this is an impure function)"
          extends Modelica.Icons.Function;
          output Integer seed "Automatically generated seed";

          external "C" seed = ModelicaRandom_automaticGlobalSeed(0.0) annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaRandom.h\"", Library="ModelicaExternalC");
         annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
seed = Utilities.<strong>automaticGlobalSeed</strong>();
</pre></blockquote>

<h4>Description</h4>
<p>Returns an automatically computed seed (Integer). Typically, this seed is computed from:</p>
<ol>
<li> The current local time by computing the number of milli-seconds up to the current hour</li>
<li> The process id (added to the first part by multiplying it with the prime number 6007).</li>
</ol>
<p>
If getTime and getPid functions are not available on the target where this Modelica function
is called, other means to compute a seed may be used.
</p>

<p>
Note, this is an impure function that returns always a different value, when it is newly called.
This function should be only called once during initialization.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>parameter</strong> Boolean useAutomaticSeed = false;
<strong>parameter</strong> Integer fixedSeed = 67867967;
<strong>final parameter</strong> Integer seed = <strong>if</strong> useAutomaticSeed <strong>then</strong>
                              Random.Utilities.automaticGlobalSeed() <strong>else</strong> fixedSeed;
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Utilities.automaticLocalSeed\">automaticLocalSeed</a>.
</p>
<h4>Note</h4>
<p>This function is impure!</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
        end automaticGlobalSeed;

        function automaticLocalSeed
          "Creates an automatic local seed from the instance name"
          extends Modelica.Icons.Function;
          input String path
            "Full path name of the instance (inquire with getInstanceName())";
          output Integer seed "Automatically generated seed";
        algorithm
          // Generate a hash value from the instance name
          seed := Modelica.Utilities.Strings.hashString(path);

         annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
seed = Utilities.<strong>automaticLocalSeed</strong>(path);
</pre></blockquote>

<h4>Description</h4>
<p>
Returns an automatically computed seed (Integer) from the hash value of
the full path name of an instance (has to be inquired in the model or block
where this function is called by the Modelica built-in operator <a href=\"https://specification.modelica.org/v3.4/Ch3.html#getinstancename\">getInstanceName()</a>).
Contrary to <a href=\"modelica://Modelica.Math.Random.Utilities.automaticGlobalSeed\">automaticGlobalSeed()</a>,
this is a pure function, that is, the same seed is returned, if an identical
path is provided.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>parameter</strong> Boolean useAutomaticLocalSeed = true;
<strong>parameter</strong> Integer fixedLocalSeed        = 10;
<strong>final parameter</strong> Integer localSeed = <strong>if</strong> useAutomaticLocalSeed <strong>then</strong>
                                   automaticLocalSeed(getInstanceName())
                                 <strong>else</strong>
                                   fixedLocalSeed;
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Utilities.automaticGlobalSeed\">automaticGlobalSeed</a>, <a href=\"modelica://Modelica.Utilities.Strings.hashString\">hashString</a> and <a href=\"https://specification.modelica.org/v3.4/Ch3.html#getinstancename\">getInstanceName</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
        end automaticLocalSeed;

        function initializeImpureRandom
          "Initializes the internal state of the impure random number generator"
          extends Modelica.Icons.Function;
          input Integer seed
            "The input seed to initialize the impure random number generator";
          output Integer id
            "Identification number to be passed as input to function impureRandom, in order that sorting is correct";
        protected
          constant Integer localSeed = 715827883
            "Since there is no local seed, a large prime number is used";
          Integer rngState[33]
            "The internal state vector of the impure random number generator";

          impure function setInternalState
            "Stores the given state vector in an external static variable"
            extends Modelica.Icons.Function;
            input Integer[33] rngState "The initial state";
            input Integer id;
            external "C" ModelicaRandom_setInternalState_xorshift1024star(rngState, size(rngState,1), id)
              annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaRandom.h\"", Library="ModelicaExternalC");
          end setInternalState;

        algorithm
          // Determine the internal state (several iterations with a generator that quickly generates good numbers
          rngState := initialStateWithXorshift64star(localSeed, seed, size(rngState, 1));
          id :=localSeed;

          // Copy the internal state into the internal C static memory
          setInternalState(rngState, id);
          annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
id = <strong>initializeImpureRandom</strong>(seed;
</pre></blockquote>

<h4>Description</h4>

<p>
Generates a hidden initial state vector for the
<a href=\"modelica://Modelica.Math.Random.Generators.Xorshift1024star\">Xorshift1024star</a>
random number generator (= xorshift1024* algorithm), from Integer input argument seed. Argument seed
can be given any value (including zero or negative number). The function returns the
dummy Integer number id. This number needs to be passed as input to function
<a href=\"modelica://Modelica.Math.Random.Utilities.impureRandom\">impureRandom</a>,
in order that the sorting order is correct (so that impureRandom is always called
after initializeImpureRandom). The function stores a reasonable initial state vector
in a C-static memory by using the
<a href=\"modelica://Modelica.Math.Random.Generators.Xorshift64star\">Xorshift64star</a>
random number generator to fill the internal state vector with 64 bit random numbers.
</p>

<h4>Example</h4>
<blockquote><pre>
  <strong>parameter</strong> Integer seed;
  Real r;
  <strong>function</strong> random = impureRandom (<strong>final</strong> id=id);
<strong>protected </strong>
  Integer id = initializeImpureRandom(seed);
<strong>equation</strong>
  // Use the random number generator
  <strong>when</strong> sample(0,0.001) <strong>then</strong>
     r = random();
  <strong>end when</strong>;
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Random.Utilities.impureRandom\">Utilities.impureRandom</a>,
<a href=\"modelica://Modelica.Math.Random.Generators\">Random.Generators</a>
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
        end initializeImpureRandom;
      annotation (Documentation(info="<html>
<p>
This package contains utility functions for the random number generators,
that are usually of no interest for the user
(they are, for example, used in package <a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a>).
</p>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end Utilities;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
        Ellipse(
          extent={{-84,84},{-24,24}},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{22,62},{82,2}},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-58,6},{2,-54}},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{26,-30},{86,-90}},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid)}), Documentation(info="<html>
<p>
This package contains low level functions for the generation of random numbers.
Usually, the functions in this package are not used directly, but are utilized
as building blocks of higher level functionality.
</p>

<p>
Package <a href=\"modelica://Modelica.Math.Random.Generators\">Math.Random.Generators</a>
contains various pseudo random number generators. These generators are used in the blocks
of package <a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a> to generate
reproducible noise signals.
Package <a href=\"modelica://Modelica.Math.Random.Utilities\">Math.Random.Utilities</a>
contains utility functions for the random number generators,
that are usually of no interest for the user
(they are, for example, used to implement the blocks in
package <a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a>).
</p>
</html>",     revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
    end Random;

    package Distributions "Library of distribution functions"
       extends Modelica.Icons.Package;

      package Uniform "Library of uniform distribution functions"
        extends Modelica.Icons.Package;

        function quantile "Quantile of uniform distribution"
          extends Modelica.Math.Distributions.Interfaces.partialQuantile;
          input Real y_min=0 "Lower limit of y" annotation (Dialog);
          input Real y_max=1 "Upper limit of y" annotation (Dialog);
        algorithm
          y := u*(y_max - y_min) + y_min;
          annotation (Inline=true,Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Uniform.<strong>quantile</strong>(u, y_min=0, y_max=1);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the inverse cumulative distribution function (= quantile) according to a <strong>uniform</strong>
distribution in a band. Input argument u must be in the range:
</p>

<blockquote>
<p>
0 &le; u &le; 1
</p>
</blockquote>

<p>
The returned number y is in the range:
</p>

<blockquote>
<p>
y_min &le; y &le; y_max
</p>
</blockquote>

<p>
Plot of the function:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/Uniform.quantile.png\">
</blockquote>

<p>
For more details, see
<a href=\"http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)\">Wikipedia</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
quantile(0.5)      // = 0.5
quantile(0.5,-1,1) // = 0
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Distributions.Uniform.density\">Uniform.density</a>,
<a href=\"modelica://Modelica.Math.Distributions.Uniform.cumulative\">Uniform.cumulative</a>.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
        end quantile;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
              -100,-100},{100,100}}), graphics={Line(
            points={{-80,-60},{-40,-60},{-40,60},{40,60},{40,-60},{80,-60}})}), Documentation(info="<html>
<p>
This package provides
</p>
<ul>
<li> probability density function (= derivative of cumulative distribution function),</li>
<li> cumulative distribution function, and</li>
<li> quantile (= inverse cumulative distribution function).</li>
</ul>
<p>
of the <strong>uniform</strong> distribution. Examples:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/Uniform.density.png\">
</blockquote>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/Uniform.cumulative.png\">
</blockquote>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/Uniform.quantile.png\">
</blockquote>

<p>
For more details of this distribution see
<a href=\"http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)\">Wikipedia</a>.
</p>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end Uniform;

      package Interfaces "Library of interfaces for distribution functions"
        extends Modelica.Icons.InterfacesPackage;

        partial function partialQuantile
          "Common interface of quantile functions (= inverse cumulative distribution functions)"
          extends Modelica.Icons.Function;
          input Real u(min=0, max=1) "Random number in the range 0 <= u <= 1";
          output Real y
            "Random number u transformed according to the given distribution";
          annotation (Documentation(info="<html>
<p>
A partial function containing the common
arguments of the quantile functions.
</p>
</html>",         revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
        end partialQuantile;
      annotation (Documentation(info="<html>
<p>
This package contains partial functions that describe the
common interface arguments of the distribution and
truncated distribution functions.
</p>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end Interfaces;
    annotation (Icon(graphics={Line(
              points={{-70,-65.953},{-66.5,-65.8975},{-63,-65.7852},{-59.5,
              -65.5674},{-56,-65.1631},{-52.5,-64.4442},{-49,-63.2213},{-45.5,
              -61.2318},{-42,-58.1385},{-38.5,-53.5468},{-35,-47.0467},{-31.5,
              -38.2849},{-28,-27.0617},{-24.5,-13.4388},{-21,2.1682},{-17.5,
              18.9428},{-14,35.695},{-10.5,50.9771},{-7,63.2797},{-3.5,
              71.2739},{0,74.047},{3.5,71.2739},{7,63.2797},{10.5,50.9771},{
              14,35.695},{17.5,18.9428},{21,2.1682},{24.5,-13.4388},{28,
              -27.0617},{31.5,-38.2849},{35,-47.0467},{38.5,-53.5468},{42,
              -58.1385},{45.5,-61.2318},{49,-63.2213},{52.5,-64.4442},{56,
              -65.1631},{59.5,-65.5674},{63,-65.7852},{66.5,-65.8975},{70,
              -65.953}},
              smooth=Smooth.Bezier)}), Documentation(info="<html>
<p>
This package provides
</p>
<ul>
<li> <a href=\"http://en.wikipedia.org/wiki/Probability_density_function\">probability density functions</a>
     (= derivative of cumulative distribution function),</li>
<li> <a href=\"http://en.wikipedia.org/wiki/Cumulative_distribution_function\">cumulative distribution functions</a>,
     and</li>
<li> <a href=\"http://en.wikipedia.org/wiki/Quantile_function\">quantiles</a>
     (= inverse cumulative distribution functions).</li>
</ul>
<p>
of different distributions.
</p>

<p>
In particular also <strong>truncated distributions</strong> are provided (see below).
The main reason to introduce
truncated distributions is to make the modeling of measurement noise easier, in order to
limit the band in which the noise can occur. For example, if a sensor is used and the
sensor signal has a noise of &plusmn; 0.1 Volt (e.g. this can be determined by using a reference
value of 0 V and inspecting the measured signal), then the sensor signal will be often the input
to an Analog-Digital converter and this converter limits the signal, say to &plusmn; 5 Volt.
Typically, the user would like to model noise within the noise band (say &plusmn; 0.1 Volt),
and often uses a normal distribution. But a normal distribution is not limited and
for a small sample time and a long simulation there might be some sample time instants
where the noise values of the normal signal is outside the &plusmn; 0.1 Volt range.
For some sensor types this is completely unrealistic (e.g. an angle sensor might
measure &plusmn; 0.1 rad, but the sensor will never add, say one revolution (6.28 rad) to it.
However, the noise model with a pure normal distribution could give such a value.
If a modeler would like to guarantee (and not to hope), that the modeled noise is
always between &plusmn; 0.1 Volt, then there are two main possibilities: (a) The noise is computed
and the result is then limited to &plusmn; 0.1 Volt, or (b) the normal distribution is slightly modified,
so that it is within the band of &plusmn; 0.1 Volt. Approach (a) is a brute force method that
changes the statistical properties of the signal in an unknown way. Approach (b)
is a \"clean\" mathematical description. The blocks in package
<a href=\"modelica://Modelica.Blocks.Noise\">Blocks.Noise</a>
give the user the freedom to choose: Either compute a normal (unlimited) noise, or
a truncated normal noise (truncated distribution).
</p>

<h4>
Details of truncated distributions
</h4>

<p>
Truncated distributions are distributions that are transformed in such a way that
either the input is within a band u_min .. u_max, or the output is within
a band y_min .. y_max.
A truncated distribution is derived from a base
distribution (e.g. from the normal distribution), by truncating its
probability density function to the desired band and adding a constant
value over this band, in order that the integral over the truncated distribution
remains one. All other properties (such as cumulative distribution function) can then be determined
in a straightforward way, provided the properties of the underlying base distribution
are available.
More details can be found, for example, in
<a href=\"http://en.wikipedia.org/wiki/Truncated_distribution\">Wikipedia</a>
(the equations from the \"Truncated Distribution\" box in the right part
of this Wikipedia article are used for this package).
</p>

<p>
When using random numbers according to a given truncated distribution,
the output of the inverse cumulative distribution function (= quantile) is restricted
to the defined band.
</p>

<p>
The truncated distribution functions are derived from the underlying distribution
functions in the following way:
</p>

<blockquote><pre>
// Original distributions
    pdf = Distributions.XXX.density(u,..);
    cdf = Distributions.XXX.cumulative(u,...);
cdf_min = Distributions.XXX.cumulative(u_min,...);
cdf_max = Distributions.XXX.cumulative(u_max,...);

// Truncated distributions
</pre></blockquote>
<blockquote>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><th><strong><em>Function</em></strong></th><th><strong><em>Transformation</em></strong></th></tr>
  <tr><td>density(u,u_min,u_max,...)</td>
      <td>= <strong>if</strong> u &ge; u_min <strong>and</strong> u&le;u_max <strong>then</strong> pdf / (cdf_max - cdf_min) <strong>else</strong> 0</td>
  </tr>
  <tr><td>cumulative(u,u_min,u_max,...)</td>
      <td>= <strong>if</strong> u &le; u_min <strong>then</strong> 0
            <strong>else if</strong> u &lt; u_max <strong>then</strong>
              (cdf - cdf_min))/(cdf_max - cdf_min)
            <strong>else</strong> 1</td>
  </tr>
  <tr><td>quantile(u,u_min,u_max,...)</td>
      <td>= Distributions.XXX.quantile( cdf_min + u*(cdf_max - cdf_min), ... )</td>
  </tr>
</table>
</blockquote>
<p>
For an example of a truncated distribution, see the following
plot of the probability density function of a normal distribution
compared with its truncated distribution:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Math/Distributions/TruncatedNormal.density.png\">
</blockquote>
</html>",     revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
    end Distributions;

  package Icons "Icons for Math"
    extends Modelica.Icons.IconsPackage;

    partial function AxisCenter
      "Basic icon for mathematical function with y-axis in the center"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{0,-80},{0,68}}, color={192,192,192}),
            Polygon(
              points={{0,90},{-8,68},{8,68},{0,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              textColor={0,0,255})}),
        Documentation(info="<html>
<p>
Icon for a mathematical function, consisting of an y-axis in the middle.
It is expected, that an x-axis is added and a plot of the function.
</p>
</html>"));
    end AxisCenter;
  end Icons;

  function asin "Inverse sine (-1 <= u <= 1)"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u "Independent variable";
    output Modelica.Units.SI.Angle y "Dependent variable y=asin(u)";

  external "builtin" y = asin(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-79.2,-72.8},{-77.6,-67.5},{-73.6,-59.4},{-66.3,
                -49.8},{-53.5,-37.3},{-30.2,-19.7},{37.4,24.8},{57.5,40.8},{68.7,
                52.7},{75.2,62.2},{77.6,67.5},{80,80}}),
          Text(
            extent={{-88,78},{-16,30}},
            textColor={192,192,192},
            textString="asin")}),
      Documentation(info="<html>
<p>
This function returns y = asin(u), with -1 &le; u &le; +1:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/asin.png\">
</p>
</html>"));
  end asin;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}), graphics={Line(points={{-80,0},{-68.7,34.2},{-61.5,53.1},
              {-55.1,66.4},{-49.4,74.6},{-43.8,79.1},{-38.2,79.8},{-32.6,76.6},{
              -26.9,69.7},{-21.3,59.4},{-14.9,44.1},{-6.83,21.2},{10.1,-30.8},{17.3,
              -50.2},{23.7,-64.2},{29.3,-73.1},{35,-78.4},{40.6,-80},{46.2,-77.6},
              {51.9,-71.5},{57.5,-61.9},{63.9,-47.2},{72,-24.8},{80,0}}, color={
              0,0,0}, smooth=Smooth.Bezier)}), Documentation(info="<html>
<p>
This package contains <strong>basic mathematical functions</strong> (such as sin(..)),
as well as functions operating on
<a href=\"modelica://Modelica.Math.Vectors\">vectors</a>,
<a href=\"modelica://Modelica.Math.Matrices\">matrices</a>,
<a href=\"modelica://Modelica.Math.Nonlinear\">nonlinear functions</a>, and
<a href=\"modelica://Modelica.Math.BooleanVectors\">Boolean vectors</a>.
</p>

<h4>Main Authors</h4>
<p><a href=\"http://www.robotic.dlr.de/Martin.Otter/\"><strong>Martin Otter</strong></a>
and <strong>Marcus Baur</strong><br>
Deutsches Zentrum f&uuml;r Luft- und Raumfahrt e.V. (DLR)<br>
Institut f&uuml;r Systemdynamik und Regelungstechnik (DLR-SR)<br>
Forschungszentrum Oberpfaffenhofen<br>
D-82234 Wessling<br>
Germany<br>
email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a>
</p>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>June 22, 2019</em>
       by Thomas Beutlich: Functions tempInterpol1/tempInterpol2 moved to ObsoleteModelica4</li>
<li><em>August 24, 2016</em>
       by Christian Kral: added wrapAngle</li>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       Function tempInterpol2 added.</li>
<li><em>Oct. 24, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Icons for icon and diagram level introduced.</li>
<li><em>June 30, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>
</html>"));
  end Math;

  package Utilities
  "Library of utility functions dedicated to scripting (operating on files, streams, strings, system)"
    extends Modelica.Icons.UtilitiesPackage;

    package Strings "Operations on strings"
      extends Modelica.Icons.FunctionsPackage;

      pure function length "Return length of string"
        extends Modelica.Icons.Function;
        input String string;
        output Integer result "Number of characters of string";
      external "C" result = ModelicaStrings_length(string) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<strong>length</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
Returns the number of characters of \"string\".
</p>
</html>"));
      end length;

      pure function compare "Compare two strings lexicographically"
        extends Modelica.Icons.Function;
        input String string1;
        input String string2;
        input Boolean caseSensitive=true "= false, if case of letters is ignored";
        output Modelica.Utilities.Types.Compare result "Result of comparison";
      external "C" result = ModelicaStrings_compare(string1, string2, caseSensitive) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = Strings.<strong>compare</strong>(string1, string2);
result = Strings.<strong>compare</strong>(string1, string2, caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Compares two strings. If the optional argument caseSensitive=false,
upper case letters are treated as if they would be lower case letters.
The result of the comparison is returned as:
</p>
<blockquote><pre>
result = Modelica.Utilities.Types.Compare.Less     // string1 &lt; string2
       = Modelica.Utilities.Types.Compare.Equal    // string1 = string2
       = Modelica.Utilities.Types.Compare.Greater  // string1 &gt; string2
</pre></blockquote>
<p>
Comparison is with regards to lexicographical order,
e.g., \"a\" &lt; \"b\";
</p>
</html>"));
      end compare;

      function isEqual "Determine whether two strings are identical"
        extends Modelica.Icons.Function;
        input String string1;
        input String string2;
        input Boolean caseSensitive=true
          "= false, if lower and upper case are ignored for the comparison";
        output Boolean identical "True, if string1 is identical to string2";
      algorithm
        identical :=compare(string1, string2, caseSensitive) == Types.Compare.Equal;
        annotation (
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<strong>isEqual</strong>(string1, string2);
Strings.<strong>isEqual</strong>(string1, string2, caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Compare whether two strings are identical,
optionally ignoring case.
</p>
</html>"));
      end isEqual;

      function isEmpty
        "Return true if a string is empty (has only white space characters)"
        extends Modelica.Icons.Function;
        input String string;
        output Boolean result "True, if string is empty";
      protected
        Integer nextIndex;
        Integer len;
      algorithm
        nextIndex := Strings.Advanced.skipWhiteSpace(string);
        len := Strings.length(string);
        if len < 1 or nextIndex > len then
          result := true;
        else
          result := false;
        end if;

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<strong>isEmpty</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
Returns true if the string has no characters or if the string consists
only of white space characters. Otherwise, false is returned.
</p>

<h4>Example</h4>
<blockquote><pre>
isEmpty(\"\");       // returns true
isEmpty(\"   \");    // returns true
isEmpty(\"  abc\");  // returns false
isEmpty(\"a\");      // returns false
</pre></blockquote>
</html>"));
      end isEmpty;

      pure function hashString "Create a hash value of a string"
        extends Modelica.Icons.Function;
        input String string "The string to create a hash from";
        output Integer hash "The hash value of string";
        external "C" hash = ModelicaStrings_hashString(string)
           annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
hash = Strings.<strong>hashString</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
Returns an Integer hash value of the provided string
(the hash can be any Integer, including zero or a negative number).
</p>

<h4>Example</h4>
<blockquote><pre>
hashString(\"this is a test\")     // =  1827717433
hashString(\"Controller.noise1\")  // = -1025762750
</pre></blockquote>
</html>",       revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td> June 22, 2015 </td>
    <td>

<table border=\"0\">
<tr><td>
         <img src=\"modelica://Modelica/Resources/Images/Logos/dlr_logo.png\" alt=\"DLR logo\">
</td><td valign=\"bottom\">
         Initial version implemented by
         A. Kl&ouml;ckner, F. v.d. Linden, D. Zimmer, M. Otter.<br>
         <a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>
</td></tr></table>
</td></tr>

</table>
</html>"));
      end hashString;

      package Advanced "Advanced scanning functions"
        extends Modelica.Icons.FunctionsPackage;

        pure function skipWhiteSpace "Scan white space"
          extends Modelica.Icons.Function;
          input String string;
          input Integer startIndex(min=1)=1;
          output Integer nextIndex;
          external "C" nextIndex = ModelicaStrings_skipWhiteSpace(string, startIndex) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
          annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
nextIndex = <strong>skipWhiteSpace</strong>(string, startIndex);
</pre></blockquote>
<h4>Description</h4>
<p>
Starts scanning of \"string\" at position \"startIndex\" and
skips white space. The function returns nextIndex = index of character
of the first non white space character.
</p>
<h4>See also</h4>
<a href=\"modelica://Modelica.Utilities.Strings.Advanced\">Strings.Advanced</a>.
</html>"));
        end skipWhiteSpace;
        annotation (Documentation(info="<html>
<h4>Library content</h4>
<p>
Package <strong>Strings.Advanced</strong> contains basic scanning
functions. These functions should be <strong>not called</strong> directly, because
it is much simpler to utilize the higher level functions \"Strings.scanXXX\".
The functions of the \"Strings.Advanced\" library provide
the basic interface in order to implement the higher level
functions in package \"Strings\".
</p>
<p>
Library \"Advanced\" provides the following functions:
</p>
<blockquote><pre>
(nextIndex, realNumber)    = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanReal\">scanReal</a>        (string, startIndex, unsigned=false);
(nextIndex, integerNumber) = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanInteger\">scanInteger</a>     (string, startIndex, unsigned=false);
(nextIndex, string2)       = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanString\">scanString</a>      (string, startIndex);
(nextIndex, identifier)    = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanIdentifier\">scanIdentifier</a>  (string, startIndex);
 nextIndex                 = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.skipWhiteSpace\">skipWhiteSpace</a>  (string, startIndex);
 nextIndex                 = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.skipLineComments\">skipLineComments</a>(string, startIndex);
</pre></blockquote>
<p>
All functions perform the following actions:
</p>
<ol>
<li> Scanning starts at character position \"startIndex\" of
     \"string\" (startIndex has a default of 1).</li>
<li> First, white space is skipped, such as blanks (\" \"), tabs (\"\\t\"), or newline (\"\\n\")</li>
<li> Afterwards, the required token is scanned.</li>
<li> If successful, on return nextIndex = index of character
     directly after the found token and the token value is returned
     as second output argument.<br>
     If not successful, on return nextIndex = startIndex.
     </li>
</ol>
<p>
The following additional rules apply for the scanning:
</p>
<ul>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanReal\">scanReal</a>:<br>
     Scans a full number including one optional leading \"+\" or \"-\" (if unsigned=false)
     according to the Modelica grammar. For example, \"+1.23e-5\", \"0.123\" are
     Real numbers, but \".1\" is not.
     Note, an Integer number, such as \"123\" is also treated as a Real number.<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanInteger\">scanInteger</a>:<br>
     Scans an Integer number including one optional leading \"+\"
     or \"-\" (if unsigned=false) according to the Modelica (and C/C++) grammar.
     For example, \"+123\", \"20\" are Integer numbers.
     Note, a Real number, such as \"123.4\" is not an Integer and
     scanInteger returns nextIndex = startIndex.<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanString\">scanString</a>:<br>
     Scans a String according to the Modelica (and C/C++) grammar, e.g.,
     \"This is a \"string\"\" is a valid string token.<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanIdentifier\">scanIdentifier</a>:<br>
     Scans a Modelica identifier, i.e., the identifier starts either
     with a letter, followed by letters, digits or \"_\".
     For example, \"w_rel\", \"T12\".<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.skipLineComments\">skipLineComments</a><br>
     Skips white space and Modelica (C/C++) line comments iteratively.
     A line comment starts with \"//\" and ends either with an
     end-of-line (\"\\n\") or the end of the \"string\".</li>
</ul>
</html>"));
      end Advanced;
      annotation (
        Documentation(info="<html>
<h4>Library content</h4>
<p>
Package <strong>Strings</strong> contains functions to manipulate strings.
</p>
<p>
In the table below an example
call to every function is given using the <strong>default</strong> options.
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><th><strong><em>Function</em></strong></th><th><strong><em>Description</em></strong></th></tr>
  <tr><td>len = <a href=\"modelica://Modelica.Utilities.Strings.length\">length</a>(string)</td>
      <td>Returns length of string</td></tr>
  <tr><td>string2 = <a href=\"modelica://Modelica.Utilities.Strings.substring\">substring</a>(string1,startIndex,endIndex)
       </td>
      <td>Returns a substring defined by start and end index</td></tr>
  <tr><td>result = <a href=\"modelica://Modelica.Utilities.Strings.repeat\">repeat</a>(n)<br>
 result = <a href=\"modelica://Modelica.Utilities.Strings.repeat\">repeat</a>(n,string)</td>
      <td>Repeat a blank or a string n times.</td></tr>
  <tr><td>result = <a href=\"modelica://Modelica.Utilities.Strings.compare\">compare</a>(string1, string2)</td>
      <td>Compares two substrings with regards to alphabetical order</td></tr>
  <tr><td>identical =
<a href=\"modelica://Modelica.Utilities.Strings.isEqual\">isEqual</a>(string1,string2)</td>
      <td>Determine whether two strings are identical</td></tr>
  <tr><td>result = <a href=\"modelica://Modelica.Utilities.Strings.count\">count</a>(string,searchString)</td>
      <td>Count the number of occurrences of a string</td></tr>
  <tr>
<td>index = <a href=\"modelica://Modelica.Utilities.Strings.find\">find</a>(string,searchString)</td>
      <td>Find first occurrence of a string in another string</td></tr>
<tr>
<td>index = <a href=\"modelica://Modelica.Utilities.Strings.findLast\">findLast</a>(string,searchString)</td>
      <td>Find last occurrence of a string in another string</td></tr>
  <tr><td>string2 = <a href=\"modelica://Modelica.Utilities.Strings.replace\">replace</a>(string,searchString,replaceString)</td>
      <td>Replace one or all occurrences of a string</td></tr>
  <tr><td>stringVector2 = <a href=\"modelica://Modelica.Utilities.Strings.sort\">sort</a>(stringVector1)</td>
      <td>Sort vector of strings in alphabetic order</td></tr>
  <tr><td>hash = <a href=\"modelica://Modelica.Utilities.Strings.hashString\">hashString</a>(string)</td>
      <td>Create a hash value of a string</td></tr>
  <tr><td>(token, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanToken\">scanToken</a>(string,startIndex)</td>
      <td>Scan for a token (Real/Integer/Boolean/String/Identifier/Delimiter/NoToken)</td></tr>
  <tr><td>(number, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanReal\">scanReal</a>(string,startIndex)</td>
      <td>Scan for a Real constant</td></tr>
  <tr><td>(number, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanInteger\">scanInteger</a>(string,startIndex)</td>
      <td>Scan for an Integer constant</td></tr>
  <tr><td>(boolean, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanBoolean\">scanBoolean</a>(string,startIndex)</td>
      <td>Scan for a Boolean constant</td></tr>
  <tr><td>(string2, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanString\">scanString</a>(string,startIndex)</td>
      <td>Scan for a String constant</td></tr>
  <tr><td>(identifier, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanIdentifier\">scanIdentifier</a>(string,startIndex)</td>
      <td>Scan for an identifier</td></tr>
  <tr><td>(delimiter, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanDelimiter\">scanDelimiter</a>(string,startIndex)</td>
      <td>Scan for delimiters</td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Strings.scanNoToken\">scanNoToken</a>(string,startIndex)</td>
      <td>Check that remaining part of string consists solely of<br>
          white space or line comments (\"// ...\\n\").</td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Strings.syntaxError\">syntaxError</a>(string,index,message)</td>
      <td> Print a \"syntax error message\" as well as a string and the<br>
           index at which scanning detected an error</td></tr>
</table>
<p>
The functions \"compare\", \"isEqual\", \"count\", \"find\", \"findLast\", \"replace\", \"sort\"
have the optional
input argument <strong>caseSensitive</strong> with default <strong>true</strong>.
If <strong>false</strong>, the operation is carried out without taking
into account whether a character is upper or lower case.
</p>
</html>"));
    end Strings;

    package Types "Type definitions used in package Modelica.Utilities"
      extends Modelica.Icons.TypesPackage;

      type Compare = enumeration(
          Less "String 1 is lexicographically less than string 2",
          Equal "String 1 is identical to string 2",
          Greater "String 1 is lexicographically greater than string 2")
        "Enumeration defining comparison of two strings";
      annotation (Documentation(info="<html>
<p>
This package contains type definitions used in Modelica.Utilities.
</p>

</html>"));
    end Types;
      annotation (
  Documentation(info="<html>
<p>
This package contains Modelica <strong>functions</strong> that are
especially suited for <strong>scripting</strong>. The functions might
be used to work with strings, read data from file, write data
to file or copy, move and remove files.
</p>
<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Utilities.UsersGuide\">Modelica.Utilities.User's Guide</a>
     discusses the most important aspects of this library.</li>
<li> <a href=\"modelica://Modelica.Utilities.Examples\">Modelica.Utilities.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
The following main sublibraries are available:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Utilities.Files\">Files</a>
     provides functions to operate on files and directories, e.g.,
     to copy, move, remove files.</li>
<li> <a href=\"modelica://Modelica.Utilities.Streams\">Streams</a>
     provides functions to read from files and write to files.</li>
<li> <a href=\"modelica://Modelica.Utilities.Strings\">Strings</a>
     provides functions to operate on strings. E.g.
     substring, find, replace, sort, scanToken.</li>
<li> <a href=\"modelica://Modelica.Utilities.System\">System</a>
     provides functions to interact with the environment.
     E.g., get or set the working directory or environment
     variables and to send a command to the default shell.</li>
</ul>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>"));
  end Utilities;

  package Constants
  "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;
    import Modelica.Units.NonSI;

    final constant Real pi=2*Modelica.Math.asin(1.0);

    final constant Real eps=ModelicaServices.Machine.eps
      "Biggest number such that 1.0 + eps = 1.0";

    final constant Real small=ModelicaServices.Machine.small
      "Smallest number such that small and -small are representable on the machine";

    final constant Real inf=ModelicaServices.Machine.inf
      "Biggest Real number such that inf and -inf are representable on the machine";

    final constant SI.Acceleration g_n=9.80665
      "Standard acceleration of gravity on earth";

    final constant SI.ElectricCharge q = 1.602176634e-19 "Elementary charge";

    final constant SI.FaradayConstant F = q*N_A
      "Faraday constant, C/mol";

    final constant Real k(final unit="J/K") = 1.380649e-23
      "Boltzmann constant";

    final constant Real R(final unit="J/(mol.K)") = k*N_A
      "Molar gas constant";

    final constant Real N_A(final unit="1/mol") = 6.02214076e23
      "Avogadro constant";

    final constant NonSI.Temperature_degC T_zero=-273.15
      "Absolute zero temperature";
    annotation (
      Documentation(info="<html>
<p>
This package provides often needed constants from mathematics, machine
dependent constants and constants from nature. The latter constants
(name, value, description) are from the following source (based on the second source):
</p>
<dl>
<dt>Michael Stock, Richard Davis, Estefan&iacute;a de Mirand&eacute;s and Martin J T Milton:</dt>
<dd><strong>The revision of the SI-the result of three decades of progress in metrology</strong> in Metrologia, Volume 56, Number 2.
<a href= \"https://iopscience.iop.org/article/10.1088/1681-7575/ab0013/pdf\">https://iopscience.iop.org/article/10.1088/1681-7575/ab0013/pdf</a>, 2019.
</dd>
</dl>
<dl>
<dt>D B Newell, F Cabiati, J Fischer, K Fujii, S G Karshenboim, H S Margolis , E de Mirand&eacute;s, P J Mohr, F Nez, K Pachucki, T J Quinn, B N Taylor, M Wang, B M Wood and Z Zhang:</dt>
<dd><strong>The CODATA 2017 values of h, e, k, and NA for the revision of the SI</strong> in Metrologia, Volume 55, Number 1.
<a href= \"https://iopscience.iop.org/article/10.1088/1681-7575/aa950a/pdf\">https://iopscience.iop.org/article/10.1088/1681-7575/aa950a/pdf</a>, 2017.
</dd>
</dl>
<p>BIPM is Bureau International des Poids et Mesures (they publish the SI-standard).</p>
<p>CODATA is the Committee on Data for Science and Technology.</p>

<dl>
<dt><strong>Main Author:</strong></dt>
<dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
    Deutsches Zentrum f&uuml;r Luft und Raumfahrt e. V. (DLR)<br>
    Oberpfaffenhofen<br>
    Postfach 1116<br>
    D-82230 We&szlig;ling<br>
    email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
</dl>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>Dec 4, 2019</em>
       by Thomas Beutlich:<br>
       Constant G updated according to 2018 CODATA value.</li>
<li><em>Mar 25, 2019</em>
       by Hans Olsson:<br>
       Constants updated according to 2017 CODATA values and new SI-standard.</li>
<li><em>Nov 4, 2015</em>
       by Thomas Beutlich:<br>
       Constants updated according to 2014 CODATA values.</li>
<li><em>Nov 8, 2004</em>
       by Christian Schweiger:<br>
       Constants updated according to 2002 CODATA values.</li>
<li><em>Dec 9, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Constants updated according to 1998 CODATA values. Using names, values
       and description text from this source. Included magnetic and
       electric constant.</li>
<li><em>Sep 18, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Constants eps, inf, small introduced.</li>
<li><em>Nov 15, 1997</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>
</html>"),
      Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
        Polygon(
          origin={-9.2597,25.6673},
          fillColor={102,102,102},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{48.017,11.336},{48.017,11.336},{10.766,11.336},{-25.684,10.95},{-34.944,-15.111},{-34.944,-15.111},{-32.298,-15.244},{-32.298,-15.244},{-22.112,0.168},{11.292,0.234},{48.267,-0.097},{48.267,-0.097}},
          smooth=Smooth.Bezier),
        Polygon(
          origin={-19.9923,-8.3993},
          fillColor={102,102,102},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{3.239,37.343},{3.305,37.343},{-0.399,2.683},{-16.936,-20.071},{-7.808,-28.604},{6.811,-22.519},{9.986,37.145},{9.986,37.145}},
          smooth=Smooth.Bezier),
        Polygon(
          origin={23.753,-11.5422},
          fillColor={102,102,102},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{-10.873,41.478},{-10.873,41.478},{-14.048,-4.162},{-9.352,-24.8},{7.912,-24.469},{16.247,0.27},{16.247,0.27},{13.336,0.071},{13.336,0.071},{7.515,-9.983},{-3.134,-7.271},{-2.671,41.214},{-2.671,41.214}},
          smooth=Smooth.Bezier)}));
  end Constants;

  package Icons "Library of icons"
    extends Icons.Package;

    partial package ExamplesPackage
    "Icon for packages containing runnable examples"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(
              origin={8.0,14.0},
              lineColor={78,138,73},
              fillColor={78,138,73},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-58.0,46.0},{42.0,-14.0},{-58.0,-74.0},{-58.0,46.0}})}), Documentation(info="<html>
<p>This icon indicates a package that contains executable examples.</p>
</html>"));
    end ExamplesPackage;

    partial package Package "Icon for standard packages"
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              lineColor={200,200,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0),
            Rectangle(
              lineColor={128,128,128},
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0)}), Documentation(info="<html>
<p>Standard package icon.</p>
</html>"));
    end Package;

    partial package InterfacesPackage "Icon for packages containing interfaces"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(origin={20.0,0.0},
              lineColor={64,64,64},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              points={{-10.0,70.0},{10.0,70.0},{40.0,20.0},{80.0,20.0},{80.0,-20.0},{40.0,-20.0},{10.0,-70.0},{-10.0,-70.0}}),
            Polygon(fillColor={102,102,102},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-100.0,20.0},{-60.0,20.0},{-30.0,70.0},{-10.0,70.0},{-10.0,-70.0},{-30.0,-70.0},{-60.0,-20.0},{-100.0,-20.0}})}),
                                Documentation(info="<html>
<p>This icon indicates packages containing interfaces.</p>
</html>"));
    end InterfacesPackage;

    partial package SourcesPackage "Icon for packages containing sources"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(origin={23.3333,0.0},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-23.333,30.0},{46.667,0.0},{-23.333,-30.0}}),
            Rectangle(
              fillColor = {128,128,128},
              pattern = LinePattern.None,
              fillPattern = FillPattern.Solid,
              extent = {{-70,-4.5},{0,4.5}})}),
                                Documentation(info="<html>
<p>This icon indicates a package which contains sources.</p>
</html>"));
    end SourcesPackage;

    partial package SensorsPackage "Icon for packages containing sensors"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(origin={0.0,-30.0},
              fillColor={255,255,255},
              extent={{-90.0,-90.0},{90.0,90.0}},
              startAngle=20.0,
              endAngle=160.0),
            Ellipse(origin={0.0,-30.0},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-20.0,-20.0},{20.0,20.0}}),
            Line(origin={0.0,-30.0},
              points={{0.0,60.0},{0.0,90.0}}),
            Ellipse(origin={-0.0,-30.0},
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-10.0,-10.0},{10.0,10.0}}),
            Polygon(
              origin={-0.0,-30.0},
              rotation=-35.0,
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-7.0,0.0},{-3.0,85.0},{0.0,90.0},{3.0,85.0},{7.0,0.0}})}),
                                Documentation(info="<html>
<p>This icon indicates a package containing sensors.</p>
</html>"));
    end SensorsPackage;

    partial package UtilitiesPackage "Icon for utility packages"
      extends Modelica.Icons.Package;
       annotation (Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
      Polygon(
        origin={1.3835,-4.1418},
        rotation=45.0,
        fillColor={64,64,64},
        pattern=LinePattern.None,
        fillPattern=FillPattern.Solid,
        points={{-15.0,93.333},{-15.0,68.333},{0.0,58.333},{15.0,68.333},{15.0,93.333},{20.0,93.333},{25.0,83.333},{25.0,58.333},{10.0,43.333},{10.0,-41.667},{25.0,-56.667},{25.0,-76.667},{10.0,-91.667},{0.0,-91.667},{0.0,-81.667},{5.0,-81.667},{15.0,-71.667},{15.0,-61.667},{5.0,-51.667},{-5.0,-51.667},{-15.0,-61.667},{-15.0,-71.667},{-5.0,-81.667},{0.0,-81.667},{0.0,-91.667},{-10.0,-91.667},{-25.0,-76.667},{-25.0,-56.667},{-10.0,-41.667},{-10.0,43.333},{-25.0,58.333},{-25.0,83.333},{-20.0,93.333}}),
      Polygon(
        origin={10.1018,5.218},
        rotation=-45.0,
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid,
        points={{-15.0,87.273},{15.0,87.273},{20.0,82.273},{20.0,27.273},{10.0,17.273},{10.0,7.273},{20.0,2.273},{20.0,-2.727},{5.0,-2.727},{5.0,-77.727},{10.0,-87.727},{5.0,-112.727},{-5.0,-112.727},{-10.0,-87.727},{-5.0,-77.727},{-5.0,-2.727},{-20.0,-2.727},{-20.0,2.273},{-10.0,7.273},{-10.0,17.273},{-20.0,27.273},{-20.0,82.273}})}),
      Documentation(info="<html>
<p>This icon indicates a package containing utility classes.</p>
</html>"));
    end UtilitiesPackage;

    partial package TypesPackage
    "Icon for packages containing type definitions"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-12.167,-23},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{12.167,65},{14.167,93},{36.167,89},{24.167,20},{4.167,-30},
                  {14.167,-30},{24.167,-30},{24.167,-40},{-5.833,-50},{-15.833,
                  -30},{4.167,20},{12.167,65}},
              smooth=Smooth.Bezier), Polygon(
              origin={2.7403,1.6673},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{49.2597,22.3327},{31.2597,24.3327},{7.2597,18.3327},{-26.7403,
                10.3327},{-46.7403,14.3327},{-48.7403,6.3327},{-32.7403,0.3327},{-6.7403,
                4.3327},{33.2597,14.3327},{49.2597,14.3327},{49.2597,22.3327}},
              smooth=Smooth.Bezier)}));
    end TypesPackage;

    partial package FunctionsPackage "Icon for packages containing functions"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
              Text(
                textColor={128,128,128},
                extent={{-90,-90},{90,90}},
                textString="f")}));
    end FunctionsPackage;

    partial package IconsPackage "Icon for packages containing icons"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}})}));
    end IconsPackage;

    partial package InternalPackage
    "Icon for an internal package (indicating that the package should not be directly utilized by user)"
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            lineColor={215,215,215},
            fillColor={255,255,255},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25),
          Rectangle(
            lineColor={215,215,215},
            extent={{-100,-100},{100,100}},
            radius=25),
          Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={215,215,215},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-55,55},{55,-55}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-60,14},{60,-14}},
            lineColor={215,215,215},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            rotation=45)}),
      Documentation(info="<html>

<p>
This icon shall be used for a package that contains internal classes not to be
directly utilized by a user.
</p>
</html>"));
    end InternalPackage;

    partial package MaterialPropertiesPackage
    "Icon for package containing property classes"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              lineColor={102,102,102},
              fillColor={204,204,204},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Sphere,
              extent={{-60.0,-60.0},{60.0,60.0}})}),
                                Documentation(info="<html>
<p>This icon indicates a package that contains properties</p>
</html>"));
    end MaterialPropertiesPackage;

    partial class RoundSensor "Icon representing a round measurement device"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Ellipse(
              fillColor={245,245,245},
              fillPattern=FillPattern.Solid,
              extent={{-70.0,-70.0},{70.0,70.0}}),
            Line(points={{0.0,70.0},{0.0,40.0}}),
            Line(points={{22.9,32.8},{40.2,57.3}}),
            Line(points={{-22.9,32.8},{-40.2,57.3}}),
            Line(points={{37.6,13.7},{65.8,23.9}}),
            Line(points={{-37.6,13.7},{-65.8,23.9}}),
            Ellipse(
              lineColor={64,64,64},
              fillColor={255,255,255},
              extent={{-12.0,-12.0},{12.0,12.0}}),
            Polygon(
              rotation=-17.5,
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-5.0,0.0},{-2.0,60.0},{0.0,65.0},{2.0,60.0},{5.0,0.0}}),
            Ellipse(
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-7.0,-7.0},{7.0,7.0}})}),
        Documentation(info="<html>
<p>
This icon is designed for a <strong>rotational sensor</strong> model.
</p>
</html>"));
    end RoundSensor;

    partial function Function "Icon for functions"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Text(
              textColor={0,0,255},
              extent={{-150,105},{150,145}},
              textString="%name"),
            Ellipse(
              lineColor = {108,88,49},
              fillColor = {255,215,136},
              fillPattern = FillPattern.Solid,
              extent = {{-100,-100},{100,100}}),
            Text(
              textColor={108,88,49},
              extent={{-90.0,-90.0},{90.0,90.0}},
              textString="f")}),
    Documentation(info="<html>
<p>This icon indicates Modelica functions.</p>
</html>"));
    end Function;

    partial record Record "Icon for records"

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
            Text(
              textColor={0,0,255},
              extent={{-150,60},{150,100}},
              textString="%name"),
            Rectangle(
              origin={0.0,-25.0},
              lineColor={64,64,64},
              fillColor={255,215,136},
              fillPattern=FillPattern.Solid,
              extent={{-100.0,-75.0},{100.0,75.0}},
              radius=25.0),
            Line(
              points={{-100.0,0.0},{100.0,0.0}},
              color={64,64,64}),
            Line(
              origin={0.0,-50.0},
              points={{-100.0,0.0},{100.0,0.0}},
              color={64,64,64}),
            Line(
              origin={0.0,-25.0},
              points={{0.0,75.0},{0.0,-75.0}},
              color={64,64,64})}), Documentation(info="<html>
<p>
This icon is indicates a record.
</p>
</html>"));
    end Record;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}})}), Documentation(info="<html>
<p>This package contains definitions for the graphical layout of components which may be used in different libraries. The icons can be utilized by inheriting them in the desired class using &quot;extends&quot; or by directly copying the &quot;icon&quot; layer.</p>

<h4>Main Authors</h4>

<dl>
<dt><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a></dt>
    <dd>Deutsches Zentrum fuer Luft und Raumfahrt e.V. (DLR)</dd>
    <dd>Oberpfaffenhofen</dd>
    <dd>Postfach 1116</dd>
    <dd>D-82230 Wessling</dd>
    <dd>email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
<dt>Christian Kral</dt>

    <dd>  <a href=\"https://christiankral.net/\">Electric Machines, Drives and Systems</a><br>
</dd>
    <dd>1060 Vienna, Austria</dd>
    <dd>email: <a href=\"mailto:dr.christian.kral@gmail.com\">dr.christian.kral@gmail.com</a></dd>
<dt>Johan Andreasson</dt>
    <dd><a href=\"https://www.modelon.com/\">Modelon AB</a></dd>
    <dd>Ideon Science Park</dd>
    <dd>22370 Lund, Sweden</dd>
    <dd>email: <a href=\"mailto:johan.andreasson@modelon.se\">johan.andreasson@modelon.se</a></dd>
</dl>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>"));
  end Icons;

  package Units "Library of type and unit definitions"
    extends Modelica.Icons.Package;

    package SI "Library of SI unit definitions"
      extends Modelica.Icons.Package;

      type Angle = Real (
          final quantity="Angle",
          final unit="rad",
          displayUnit="deg");

      type Volume = Real (final quantity="Volume", final unit="m3");

      type Time = Real (final quantity="Time", final unit="s");

      type Velocity = Real (final quantity="Velocity", final unit="m/s");

      type Acceleration = Real (final quantity="Acceleration", final unit="m/s2");

      type Period = Real (final quantity="Time", final unit="s");

      type Frequency = Real (final quantity="Frequency", final unit="Hz");

      type Mass = Real (
          quantity="Mass",
          final unit="kg",
          min=0);

      type Density = Real (
          final quantity="Density",
          final unit="kg/m3",
          displayUnit="g/cm3",
          min=0.0);

      type SpecificVolume = Real (
          final quantity="SpecificVolume",
          final unit="m3/kg",
          min=0.0);

      type Pressure = Real (
          final quantity="Pressure",
          final unit="Pa",
          displayUnit="bar");

      type AbsolutePressure = Pressure (min=0.0, nominal = 1e5);

      type DynamicViscosity = Real (
          final quantity="DynamicViscosity",
          final unit="Pa.s",
          min=0);

      type Power = Real (final quantity="Power", final unit="W");

      type MassFlowRate = Real (quantity="MassFlowRate", final unit="kg/s");

      type VolumeFlowRate = Real (final quantity="VolumeFlowRate", final unit=
              "m3/s");

      type ThermodynamicTemperature = Real (
          final quantity="ThermodynamicTemperature",
          final unit="K",
          min = 0.0,
          start = 288.15,
          nominal = 300,
          displayUnit="degC")
        "Absolute temperature (use type TemperatureDifference for relative temperatures)" annotation(absoluteValue=true);

      type Temperature = ThermodynamicTemperature;

      type Compressibility = Real (final quantity="Compressibility", final unit=
              "1/Pa");

      type IsothermalCompressibility = Compressibility;

      type Heat = Real (final quantity="Energy", final unit="J");

      type HeatFlowRate = Real (final quantity="Power", final unit="W");

      type ThermalConductivity = Real (final quantity="ThermalConductivity", final unit=
                 "W/(m.K)");

      type SpecificHeatCapacity = Real (final quantity="SpecificHeatCapacity",
            final unit="J/(kg.K)");

      type RatioOfSpecificHeatCapacities = Real (final quantity=
              "RatioOfSpecificHeatCapacities", final unit="1");

      type SpecificEntropy = Real (final quantity="SpecificEntropy",
                                   final unit="J/(kg.K)");

      type Enthalpy = Heat;

      type SpecificEnergy = Real (final quantity="SpecificEnergy",
                                  final unit="J/kg");

      type SpecificEnthalpy = SpecificEnergy;

      type DerDensityByEnthalpy = Real (final unit="kg.s2/m5");

      type DerDensityByPressure = Real (final unit="s2/m2");

      type DerDensityByTemperature = Real (final unit="kg/(m3.K)");

      type ElectricCurrent = Real (final quantity="ElectricCurrent", final unit="A");

      type Current = ElectricCurrent;

      type ElectricCharge = Real (final quantity="ElectricCharge", final unit="C");

      type ElectricPotential = Real (final quantity="ElectricPotential", final unit=
             "V");

      type MolarMass = Real (final quantity="MolarMass", final unit="kg/mol", min=0);

      type MolarVolume = Real (final quantity="MolarVolume", final unit="m3/mol", min=0);

      type MolarEnergy = Real (final quantity="MolarEnergy", final unit="J/mol", nominal=2e4);

      type MolarHeatCapacity = Real (final quantity="MolarHeatCapacity", final unit="J/(mol.K)");

      type MolarEntropy = Real (final quantity="MolarEntropy", final unit="J/(mol.K)");

      type MolarEnthalpy = MolarEnergy;

      type MolarFlowRate = Real (final quantity="MolarFlowRate", final unit="mol/s");

      type MassConcentration = Real (final quantity="MassConcentration", final unit=
             "kg/m3");

      type MassFraction = Real (final quantity="MassFraction", final unit="1",
                                min=0, max=1);

      type Concentration = Real (final quantity="Concentration", final unit=
              "mol/m3");

      type MoleFraction = Real (final quantity="MoleFraction", final unit="1",
                                min = 0, max = 1);

      type ChemicalPotential = Real (final quantity="ChemicalPotential", final unit=
             "J/mol");

      type ActivityCoefficient = Real (final quantity="ActivityCoefficient", final unit=
                 "1");

      type ChargeNumberOfIon = Real (final quantity="ChargeNumberOfIon", final unit=
             "1");

      type FaradayConstant = Real (final quantity="FaradayConstant", final unit=
              "C/mol");

      type PrandtlNumber = Real (final quantity="PrandtlNumber", final unit="1");
      annotation (Icon(graphics={Text(
              extent={{-80,80},{80,-78}},
              textColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.None,
              fontName="serif",
              textString="SI",
              textStyle={TextStyle.Italic})}),
                                       Documentation(info="<html>
<p>This package provides predefined types based on the international standard
on units.
</p>
<p>
For an introduction to the conventions used in this package, have a look at:
<a href=\"modelica://Modelica.Units.UsersGuide.Conventions\">Conventions</a>.
</p>
</html>"));
    end SI;

    package NonSI "Type definitions of non SI and other units"
      extends Modelica.Icons.Package;

      type Temperature_degC = Real (final quantity="ThermodynamicTemperature",
            final unit="degC")
        "Absolute temperature in degree Celsius (for relative temperature use Modelica.Units.SI.TemperatureDifference)" annotation(absoluteValue=true);

      type Pressure_bar = Real (final quantity="Pressure", final unit="bar")
        "Absolute pressure in bar";
      annotation (Documentation(info="<html>
<p>
This package provides predefined types, such as <strong>Angle_deg</strong> (angle in
degree), <strong>AngularVelocity_rpm</strong> (angular velocity in revolutions per
minute) or <strong>Temperature_degF</strong> (temperature in degree Fahrenheit),
which are in common use but are not part of the international standard on
units according to ISO 31-1992 \"General principles concerning quantities,
units and symbols\" and ISO 1000-1992 \"SI units and recommendations for
the use of their multiples and of certain other units\".</p>
<p>If possible, the types in this package should not be used. Use instead
types of package <code>Modelica.Units.SI</code>. For more information on units, see also
the book of Francois Cardarelli <strong>Scientific Unit Conversion - A
Practical Guide to Metrication</strong> (Springer 1997).</p>
</html>"), Icon(coordinateSystem(extent={{-100,-100},{100,100}}), graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{-60,10},{-40,-10}},
              lineColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{40,10},{60,-10}},
              lineColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.Solid)}));
    end NonSI;

    package Conversions
    "Conversion functions to/from non SI units and type definitions of non SI units"
      extends Modelica.Icons.Package;

      function to_degC "Convert from kelvin to degree Celsius"
        extends Modelica.Units.Icons.Conversion;
        input SI.Temperature Kelvin "Value in kelvin";
        output Modelica.Units.NonSI.Temperature_degC Celsius "Value in degree Celsius";
      algorithm
        Celsius := Kelvin + Modelica.Constants.T_zero;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-20,100},{-100,20}},
                textString="K"), Text(
                extent={{100,-20},{20,-100}},
                textString="degC")}));
      end to_degC;

      function from_degC "Convert from degree Celsius to kelvin"
        extends Modelica.Units.Icons.Conversion;
        input Modelica.Units.NonSI.Temperature_degC Celsius "Value in degree Celsius";
        output SI.Temperature Kelvin "Value in kelvin";
      algorithm
        Kelvin := Celsius - Modelica.Constants.T_zero;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-20,100},{-100,20}},
                textString="degC"), Text(
                extent={{100,-20},{20,-100}},
                textString="K")}));
      end from_degC;

      function to_bar "Convert from Pascal to bar"
        extends Modelica.Units.Icons.Conversion;
        input SI.Pressure Pa "Value in Pascal";
        output Modelica.Units.NonSI.Pressure_bar bar "Value in bar";
      algorithm
        bar := Pa/1e5;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-12,100},{-100,56}},
                textString="Pa"), Text(
                extent={{98,-52},{-4,-100}},
                textString="bar")}));
      end to_bar;
      annotation (Documentation(info="<html>
<p>This package provides conversion functions from the non SI Units
defined in package <code>Modelica.Units.NonSI</code> to the
corresponding SI Units defined in package <code>Modelica.Units.SI</code> and vice
versa. It is recommended to use these functions in the following
way (note, that all functions have one Real input and one Real output
argument):</p>
<blockquote><pre>
<strong>import</strong> Modelica.Units.SI;
<strong>import</strong> Modelica.Units.Conversions.{from_degC, from_deg, from_rpm};
   ...
<strong>parameter</strong> SI.Temperature     T   = from_degC(25);   // convert 25 degree Celsius to kelvin
<strong>parameter</strong> SI.Angle           phi = from_deg(180);   // convert 180 degree to radian
<strong>parameter</strong> SI.AngularVelocity w   = from_rpm(3600);  // convert 3600 revolutions per minutes
                                                   // to radian per seconds
</pre></blockquote>

</html>"),     Icon(graphics={
            Polygon(
              points={{80,0},{20,20},{20,-20},{80,0}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,0},{20,0}}, color={191,0,0})}));
    end Conversions;

    package Icons "Icons for Units"
      extends Modelica.Icons.IconsPackage;

      partial function Conversion "Base icon for conversion functions"

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={
              Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={191,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,0},{30,0}}, color={191,0,0}),
              Polygon(
                points={{90,0},{30,20},{30,-20},{90,0}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-115,155},{115,105}},
                textString="%name",
                textColor={0,0,255})}));
      end Conversion;
    end Icons;
    annotation (Icon(graphics={
        Polygon(
          fillColor = {128,128,128},
          pattern = LinePattern.None,
          fillPattern = FillPattern.Solid,
          points = {{-80,-40},{-80,-40},{-55,50},{-52.5,62.5},{-65,60},{-65,65},{-35,77.5},{-32.5,60},{-50,0},{-50,0},{-30,15},{-20,27.5},{-32.5,27.5},{-32.5,27.5},{-32.5,32.5},{-32.5,32.5},{2.5,32.5},{2.5,32.5},{2.5,27.5},{2.5,27.5},{-7.5,27.5},{-30,7.5},{-30,7.5},{-25,-25},{-17.5,-28.75},{-10,-25},{-5,-26.25},{-5,-32.5},{-16.25,-41.25},{-31.25,-43.75},{-40,-33.75},{-45,-5},{-45,-5},{-52.5,-10},{-52.5,-10},{-60,-40},{-60,-40}},
          smooth = Smooth.Bezier),
        Polygon(
          fillColor = {128,128,128},
          pattern = LinePattern.None,
          fillPattern = FillPattern.Solid,
          points = {{87.5,30},{62.5,30},{62.5,30},{55,33.75},{36.25,35},{16.25,25},{7.5,6.25},{11.25,-7.5},{22.5,-12.5},{22.5,-12.5},{6.25,-22.5},{6.25,-35},{16.25,-38.75},{16.25,-38.75},{21.25,-41.25},{21.25,-41.25},{45,-48.75},{47.5,-61.25},{32.5,-70},{12.5,-65},{7.5,-51.25},{21.25,-41.25},{21.25,-41.25},{16.25,-38.75},{16.25,-38.75},{6.25,-41.25},{-6.25,-50},{-3.75,-68.75},{30,-76.25},{65,-62.5},{63.75,-35},{27.5,-26.25},{22.5,-20},{27.5,-15},{27.5,-15},{30,-7.5},{30,-7.5},{27.5,-2.5},{28.75,11.25},{36.25,27.5},{47.5,30},{53.75,22.5},{51.25,8.75},{45,-6.25},{35,-11.25},{30,-7.5},{30,-7.5},{27.5,-15},{27.5,-15},{43.75,-16.25},{65,-6.25},{72.5,10},{70,20},{70,20},{80,20}},
          smooth = Smooth.Bezier)}), Documentation(info="<html>
<p>This package provides predefined types, such as <em>Mass</em>,
<em>Angle</em>, <em>Time</em>, based on the international standard
on units, e.g.,
</p>

<blockquote><pre>
<strong>type</strong> Angle = Real(<strong>final</strong> quantity = \"Angle\",
                  <strong>final</strong> unit     = \"rad\",
                  displayUnit   = \"deg\");
</pre></blockquote>

<p>
Some of the types are derived SI units that are utilized in package Modelica
(such as ComplexCurrent, which is a complex number where both the real and imaginary
part have the SI unit Ampere).
</p>

<p>
Furthermore, conversion functions from non SI-units to SI-units and vice versa
are provided in subpackage
<a href=\"modelica://Modelica.Units.Conversions\">Conversions</a>.
</p>

<p>
For an introduction how units are used in the Modelica Standard Library
with package Units, have a look at:
<a href=\"modelica://Modelica.Units.UsersGuide.HowToUseUnits\">How to use Units</a>.
</p>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>May 25, 2011</em> by Stefan Wischhusen:<br>Added molar units for energy and enthalpy.</li>
<li><em>Jan. 27, 2010</em> by Christian Kral:<br>Added complex units.</li>
<li><em>Dec. 14, 2005</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Add User&#39;s Guide and removed &quot;min&quot; values for Resistance and Conductance.</li>
<li><em>October 21, 2002</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Christian Schweiger:<br>Added new package <strong>Conversions</strong>. Corrected typo <em>Wavelenght</em>.</li>
<li><em>June 6, 2000</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Introduced the following new types<br>type Temperature = ThermodynamicTemperature;<br>types DerDensityByEnthalpy, DerDensityByPressure, DerDensityByTemperature, DerEnthalpyByPressure, DerEnergyByDensity, DerEnergyByPressure<br>Attribute &quot;final&quot; removed from min and max values in order that these values can still be changed to narrow the allowed range of values.<br>Quantity=&quot;Stress&quot; removed from type &quot;Stress&quot;, in order that a type &quot;Stress&quot; can be connected to a type &quot;Pressure&quot;.</li>
<li><em>Oct. 27, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>New types due to electrical library: Transconductance, InversePotential, Damping.</li>
<li><em>Sept. 18, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Renamed from SIunit to SIunits. Subpackages expanded, i.e., the SIunits package, does no longer contain subpackages.</li>
<li><em>Aug 12, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Type &quot;Pressure&quot; renamed to &quot;AbsolutePressure&quot; and introduced a new type &quot;Pressure&quot; which does not contain a minimum of zero in order to allow convenient handling of relative pressure. Redefined BulkModulus as an alias to AbsolutePressure instead of Stress, since needed in hydraulics.</li>
<li><em>June 29, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Bug-fix: Double definition of &quot;Compressibility&quot; removed and appropriate &quot;extends Heat&quot; clause introduced in package SolidStatePhysics to incorporate ThermodynamicTemperature.</li>
<li><em>April 8, 1998</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Astrid Jaschinski:<br>Complete ISO 31 chapters realized.</li>
<li><em>Nov. 15, 1997</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Hubertus Tummescheit:<br>Some chapters realized.</li>
</ul>
</html>"));
  end Units;
annotation (
preferredView="info",
version="4.0.0",
versionDate="2020-06-04",
dateModified = "2020-06-04 11:00:00Z",
revisionId="6626538a2 2020-06-04 19:56:34 +0200",
uses(Complex(version="4.0.0"), ModelicaServices(version="4.0.0")),
conversion(
 from(version={"3.0", "3.0.1", "3.1", "3.2", "3.2.1", "3.2.2", "3.2.3"}, script="modelica://Modelica/Resources/Scripts/Conversion/ConvertModelica_from_3.2.3_to_4.0.0.mos")),
Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
  Polygon(
    origin={-6.9888,20.048},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    points={{-93.0112,10.3188},{-93.0112,10.3188},{-73.011,24.6},{-63.011,31.221},{-51.219,36.777},{-39.842,38.629},{-31.376,36.248},{-25.819,29.369},{-24.232,22.49},{-23.703,17.463},{-15.501,25.135},{-6.24,32.015},{3.02,36.777},{15.191,39.423},{27.097,37.306},{32.653,29.633},{35.035,20.108},{43.501,28.046},{54.085,35.19},{65.991,39.952},{77.897,39.688},{87.422,33.338},{91.126,21.696},{90.068,9.525},{86.099,-1.058},{79.749,-10.054},{71.283,-21.431},{62.816,-33.337},{60.964,-32.808},{70.489,-16.14},{77.368,-2.381},{81.072,10.054},{79.749,19.05},{72.605,24.342},{61.758,23.019},{49.587,14.817},{39.003,4.763},{29.214,-6.085},{21.012,-16.669},{13.339,-26.458},{5.401,-36.777},{-1.213,-46.037},{-6.24,-53.446},{-8.092,-52.387},{-0.684,-40.746},{5.401,-30.692},{12.81,-17.198},{19.424,-3.969},{23.658,7.938},{22.335,18.785},{16.514,23.283},{8.047,23.019},{-1.478,19.05},{-11.267,11.113},{-19.734,2.381},{-29.259,-8.202},{-38.519,-19.579},{-48.044,-31.221},{-56.511,-43.392},{-64.449,-55.298},{-72.386,-66.939},{-77.678,-74.612},{-79.53,-74.083},{-71.857,-61.383},{-62.861,-46.037},{-52.278,-28.046},{-44.869,-15.346},{-38.784,-2.117},{-35.344,8.731},{-36.403,19.844},{-42.488,23.813},{-52.013,22.49},{-60.744,16.933},{-68.947,10.054},{-76.884,2.646},{-93.0112,-12.1707},{-93.0112,-12.1707}},
    smooth=Smooth.Bezier),
  Ellipse(
    origin={40.8208,-37.7602},
    fillColor={161,0,4},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    extent={{-17.8562,-17.8563},{17.8563,17.8562}})}),
Documentation(info="<html>
<p>
<img src=\"modelica://Modelica/Resources/Images/Logos/Modelica_Libraries.svg\" width=\"250\">
</p>

<p>
The package <strong>Modelica&reg;</strong> is a <strong>standardized</strong> and <strong>free</strong> package
that is developed by the \"<strong>Modelica Association Project - Libraries</strong>\".</p>
<p>
Its development is coordinated with the Modelica&reg; language from the
Modelica Association, see <a href=\"https://www.Modelica.org\">https://www.Modelica.org</a>.
It is also called <strong>Modelica Standard Library</strong>.
It provides model components in many domains that are based on
standardized interface definitions. Some typical examples are shown
in the next figure:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/UsersGuide/ModelicaLibraries.png\">
</p>

<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.UsersGuide.Overview\">Overview</a>
  provides an overview of the Modelica Standard Library
  inside the <a href=\"modelica://Modelica.UsersGuide\">User's Guide</a>.</li>
<li><a href=\"modelica://Modelica.UsersGuide.ReleaseNotes\">Release Notes</a>
 summarizes the changes of new versions of this package.</li>
<li> <a href=\"modelica://Modelica.UsersGuide.Contact\">Contact</a>
  lists the contributors of the Modelica Standard Library.</li>
<li> The <strong>Examples</strong> packages in the various libraries, demonstrate
  how to use the components of the corresponding sublibrary.</li>
</ul>

<p>
This version of the Modelica Standard Library consists of
</p>
<ul>
<li><strong>1417</strong> component models and blocks,</li>
<li><strong>512</strong> example models, and</li>
<li><strong>1219</strong> functions</li>
</ul>
<p>
that are directly usable (= number of public, non-partial, non-internal and non-obsolete classes). It is fully compliant
to <a href=\"https://modelica.org/documents/ModelicaSpec34.pdf\">Modelica Specification version 3.4</a>
and it has been tested with Modelica tools from different vendors.
</p>

<p>
<strong>Licensed by the Modelica Association under the 3-Clause BSD License</strong><br>
Copyright &copy; 1998-2020, Modelica Association and <a href=\"modelica://Modelica.UsersGuide.Contact\">contributors</a>.
</p>

<p>
<em>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the 3-Clause BSD license. For license conditions (including the disclaimer of warranty) visit <a href=\"https://modelica.org/licenses/modelica-3-clause-bsd\">https://modelica.org/licenses/modelica-3-clause-bsd</a>.</em>
</p>

<p>
<strong>Modelica&reg;</strong> is a registered trademark of the Modelica Association.
</p>
</html>"));
end Modelica;

package modelECMORespiratoryVR

  package BloodGasesTransport
  "Transport of O2 and CO2 through respiration and circulation in human body"

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
      end Parts;

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
  end BloodGasesTransport;
  annotation (uses(
      Modelica(version="4.0.0"),
      Chemical(version="1.4.0"),
      Physiolibrary(version="3.0.0-beta1")),   version="1");
end modelECMORespiratoryVR;
model modelECMORespiratoryVR_BloodGasesTransport_MeursModel2011_HemodynamicsRegulatedHR
 extends modelECMORespiratoryVR.BloodGasesTransport.MeursModel2011.HemodynamicsRegulatedHR;
  annotation(experiment(StopTime=10, __Dymola_Algorithm="Dassl"),uses(modelECMORespiratoryVR(version="1")));
end modelECMORespiratoryVR_BloodGasesTransport_MeursModel2011_HemodynamicsRegulatedHR;
