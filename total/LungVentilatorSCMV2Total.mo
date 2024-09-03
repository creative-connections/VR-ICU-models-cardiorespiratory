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

package Chemical "Physical Chemistry"
 extends Modelica.Icons.Package;

  package Substances "Definitions of substances"
      extends Modelica.Icons.Package;

    record CarbonDioxide_gas "CO2(g)"
     extends Chemical.Interfaces.IdealGas.SubstanceData(
        MolarWeight=0.044,
        DfH=-393500,
        DfG=-394400,
        Cp=37.1,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end CarbonDioxide_gas;

    record Water_gas "H2O(g)"
     extends Chemical.Interfaces.IdealGas.SubstanceData(
        MolarWeight=0.018015,
        DfH=-241830,
        DfG=-228590,
        Cp=33.6,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Water_gas;

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

    record Oxygen_gas "O2(g)"
     extends Chemical.Interfaces.IdealGas.SubstanceData(
        MolarWeight=0.032,
        DfH=0,
        DfG=0,
        Cp=29.4,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Oxygen_gas;

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

    record Nitrogen_gas "N2(g)"
       extends Chemical.Interfaces.IdealGas.SubstanceData(
          MolarWeight=0.0280134,
          DfH=0,
          DfG=0,
          Cp=29.1,
          References={
              "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Type=JANAFG&Plot=on"});
        annotation (preferredView = "info");
    end Nitrogen_gas;

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

    package IdealGas "Ideal gas with constant heat capacity"
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

        parameter Modelica.Units.SI.MolarEnergy DfG_25degC_1bar(displayUnit=
           "kJ/mol") = 0 "Obsolete parameter use DfH instead"
        annotation (Dialog(tab="Obsolete"));

        parameter Modelica.Units.SI.MolarEnergy DfH_25degC(displayUnit=
            "kJ/mol") = 0 "Obsolete parameter use DfG instead"
        annotation (Dialog(tab="Obsolete"));

       parameter Boolean SelfClustering = false "Pure substance is making clusters (weak bonds between molecules)";

        parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
        "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";                                                                    //-20000
        parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
        "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

        annotation ( preferredView = "info", Documentation(revisions="<html>
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
      "Molar enthalpy of the pure substance in electroneutral solution"
     algorithm
         //Molar enthalpy:
         // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
         //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
         //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
         molarEnthalpyElectroneutral := substanceData.DfH
           +(T-298.15)*(substanceData.Cp);
     end molarEnthalpyElectroneutral;

     redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
     algorithm
       //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
       //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

         //Molar entropy:
         // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
         // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
         molarEntropyPure := (substanceData.Cp)*log(T/298.15) - Modelica.Constants.R*log(p/100000) + ((substanceData.DfH
          - substanceData.DfG)/298.15);

         //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
         //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
     end molarEntropyPure;

     redeclare function extends molarVolumePure
      "Molar volume of the pure substance"
     algorithm
         molarVolumePure := Modelica.Constants.R*T/p; //ideal gas
     end molarVolumePure;

     redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
     algorithm
         molarHeatCapacityCp := substanceData.Cp;
     end molarHeatCapacityCp;

     redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
     algorithm
         molarMass := substanceData.MolarWeight;
     end molarMassOfBaseMolecule;

     redeclare function extends temperature "Temperature of substance from its enthalpy"
     algorithm
          T := 298.15 + (h-specificEnthalpy(substanceData,298.15,p,v,I))/specificHeatCapacityCp(substanceData);
     end temperature;

     redeclare function extends solution_temperature
      "Temperature of the solution from enthalpies os substances"
     algorithm
          T := 298.15 + (h-X*specificEnthalpy(
               substanceData,
               298.15,
               p,
               v,
               I))/(X*specificHeatCapacityCp(substanceData));
     end solution_temperature;

     redeclare function extends density
          "Return density of the substance in the solution"
     algorithm
             density := substanceData.MolarWeight * p / (Modelica.Constants.R * T);
     end density;
      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end IdealGas;
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

package Modelica "Modelica Standard Library - Version 4.0.0"
extends Modelica.Icons.Package;

  package Blocks
  "Library of basic input/output control blocks (continuous, discrete, logical, table blocks)"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;

    package Continuous
    "Library of continuous control blocks with internal states"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.Package;

      block Filter
        "Continuous low pass, high pass, band pass or band stop IIR-filter of type CriticalDamping, Bessel, Butterworth or ChebyshevI"
        import Modelica.Blocks.Continuous.Internal;

        extends Modelica.Blocks.Interfaces.SISO;

        parameter Modelica.Blocks.Types.AnalogFilter analogFilter=Modelica.Blocks.Types.AnalogFilter.CriticalDamping
          "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/ChebyshevI)";
        parameter Modelica.Blocks.Types.FilterType filterType=Modelica.Blocks.Types.FilterType.LowPass
          "Type of filter (LowPass/HighPass/BandPass/BandStop)";
        parameter Integer order(min=1) = 2 "Order of filter";
        parameter SI.Frequency f_cut "Cut-off frequency";
        parameter Real gain=1.0
          "Gain (= amplitude of frequency response at zero frequency)";
        parameter Real A_ripple(unit="dB") = 0.5
          "Pass band ripple for Chebyshev filter (otherwise not used); > 0 required"
          annotation(Dialog(enable=analogFilter==Modelica.Blocks.Types.AnalogFilter.ChebyshevI));
        parameter SI.Frequency f_min=0
          "Band of band pass/stop filter is f_min (A=-3db*gain) .. f_cut (A=-3db*gain)"
          annotation(Dialog(enable=filterType == Modelica.Blocks.Types.FilterType.BandPass or
                                   filterType == Modelica.Blocks.Types.FilterType.BandStop));
        parameter Boolean normalized=true
          "= true, if amplitude at f_cut = -3db, otherwise unmodified filter";
        parameter Modelica.Blocks.Types.Init init=Modelica.Blocks.Types.Init.SteadyState
          "Type of initialization (no init/steady state/initial state/initial output)"
          annotation(Evaluate=true, Dialog(tab="Advanced"));
        final parameter Integer nx = if filterType == Modelica.Blocks.Types.FilterType.LowPass or
                                        filterType == Modelica.Blocks.Types.FilterType.HighPass then
                                        order else 2*order;
        parameter Real x_start[nx] = zeros(nx) "Initial or guess values of states"
          annotation(Dialog(tab="Advanced"));
        parameter Real y_start = 0 "Initial value of output"
          annotation(Dialog(tab="Advanced"));
        parameter Real u_nominal = 1.0
          "Nominal value of input (used for scaling the states)"
        annotation(Dialog(tab="Advanced"));
        Modelica.Blocks.Interfaces.RealOutput x[nx] "Filter states";

      protected
        parameter Integer ncr = if analogFilter == Modelica.Blocks.Types.AnalogFilter.CriticalDamping then
                                   order else mod(order,2);
        parameter Integer nc0 = if analogFilter == Modelica.Blocks.Types.AnalogFilter.CriticalDamping then
                                   0 else integer(order/2);
        parameter Integer na = if filterType == Modelica.Blocks.Types.FilterType.BandPass or
                                  filterType == Modelica.Blocks.Types.FilterType.BandStop then order else
                               if analogFilter == Modelica.Blocks.Types.AnalogFilter.CriticalDamping then
                                  0 else integer(order/2);
        parameter Integer nr = if filterType == Modelica.Blocks.Types.FilterType.BandPass or
                                  filterType == Modelica.Blocks.Types.FilterType.BandStop then 0 else
                               if analogFilter == Modelica.Blocks.Types.AnalogFilter.CriticalDamping then
                                  order else mod(order,2);

        // Coefficients of prototype base filter (low pass filter with w_cut = 1 rad/s)
        parameter Real cr[ncr](each fixed=false);
        parameter Real c0[nc0](each fixed=false);
        parameter Real c1[nc0](each fixed=false);

        // Coefficients for differential equations.
        parameter Real r[nr](each fixed=false);
        parameter Real a[na](each fixed=false);
        parameter Real b[na](each fixed=false);
        parameter Real ku[na](each fixed=false);
        parameter Real k1[if filterType == Modelica.Blocks.Types.FilterType.LowPass then 0 else na](
                       each fixed = false);
        parameter Real k2[if filterType == Modelica.Blocks.Types.FilterType.LowPass then 0 else na](
                       each fixed = false);

        // Auxiliary variables
        Real uu[na+nr+1];

      initial equation
         if analogFilter == Modelica.Blocks.Types.AnalogFilter.CriticalDamping then
            cr = Internal.Filter.base.CriticalDamping(order, normalized);
         elseif analogFilter == Modelica.Blocks.Types.AnalogFilter.Bessel then
            (cr,c0,c1) = Internal.Filter.base.Bessel(order, normalized);
         elseif analogFilter == Modelica.Blocks.Types.AnalogFilter.Butterworth then
            (cr,c0,c1) = Internal.Filter.base.Butterworth(order, normalized);
         elseif analogFilter == Modelica.Blocks.Types.AnalogFilter.ChebyshevI then
            (cr,c0,c1) = Internal.Filter.base.ChebyshevI(order, A_ripple, normalized);
         end if;

         if filterType == Modelica.Blocks.Types.FilterType.LowPass then
            (r,a,b,ku) = Internal.Filter.roots.lowPass(cr,c0,c1,f_cut);
         elseif filterType == Modelica.Blocks.Types.FilterType.HighPass then
            (r,a,b,ku,k1,k2) = Internal.Filter.roots.highPass(cr,c0,c1,f_cut);
         elseif filterType == Modelica.Blocks.Types.FilterType.BandPass then
            (a,b,ku,k1,k2) = Internal.Filter.roots.bandPass(cr,c0,c1,f_min,f_cut);
         elseif filterType == Modelica.Blocks.Types.FilterType.BandStop then
            (a,b,ku,k1,k2) = Internal.Filter.roots.bandStop(cr,c0,c1,f_min,f_cut);
         end if;

         if init == Modelica.Blocks.Types.Init.InitialState then
            x = x_start;
         elseif init == Modelica.Blocks.Types.Init.SteadyState then
            der(x) = zeros(nx);
         elseif init == Modelica.Blocks.Types.Init.InitialOutput then
            y = y_start;
            if nx > 1 then
               der(x[1:nx-1]) = zeros(nx-1);
            end if;
         end if;

      equation
         assert(u_nominal > 0, "u_nominal > 0 required");
         assert(filterType == Modelica.Blocks.Types.FilterType.LowPass or
                filterType == Modelica.Blocks.Types.FilterType.HighPass or
                f_min > 0, "f_min > 0 required for band pass and band stop filter");
         assert(A_ripple > 0, "A_ripple > 0 required");
         assert(f_cut > 0, "f_cut > 0 required");

         /* All filters have the same basic differential equations:
        Real poles:
           der(x) = r*x - r*u
        Complex conjugate poles:
           der(x1) = a*x1 - b*x2 + ku*u;
           der(x2) = b*x1 + a*x2;
   */
         uu[1] = u/u_nominal;
         for i in 1:nr loop
            der(x[i]) = r[i]*(x[i] - uu[i]);
         end for;
         for i in 1:na loop
            der(x[nr+2*i-1]) = a[i]*x[nr+2*i-1] - b[i]*x[nr+2*i] + ku[i]*uu[nr+i];
            der(x[nr+2*i])   = b[i]*x[nr+2*i-1] + a[i]*x[nr+2*i];
         end for;

         // The output equation is different for the different filter types
         if filterType == Modelica.Blocks.Types.FilterType.LowPass then
            /* Low pass filter
           Real poles             :  y = x
           Complex conjugate poles:  y = x2
      */
            for i in 1:nr loop
               uu[i+1] = x[i];
            end for;
            for i in 1:na loop
               uu[nr+i+1] = x[nr+2*i];
            end for;

         elseif filterType == Modelica.Blocks.Types.FilterType.HighPass then
            /* High pass filter
           Real poles             :  y = -x + u;
           Complex conjugate poles:  y = k1*x1 + k2*x2 + u;
      */
            for i in 1:nr loop
               uu[i+1] = -x[i] + uu[i];
            end for;
            for i in 1:na loop
               uu[nr+i+1] = k1[i]*x[nr+2*i-1] + k2[i]*x[nr+2*i] + uu[nr+i];
            end for;

         elseif filterType == Modelica.Blocks.Types.FilterType.BandPass then
            /* Band pass filter
           Complex conjugate poles:  y = k1*x1 + k2*x2;
      */
            for i in 1:na loop
               uu[nr+i+1] = k1[i]*x[nr+2*i-1] + k2[i]*x[nr+2*i];
            end for;

         elseif filterType == Modelica.Blocks.Types.FilterType.BandStop then
            /* Band pass filter
           Complex conjugate poles:  y = k1*x1 + k2*x2 + u;
      */
            for i in 1:na loop
               uu[nr+i+1] = k1[i]*x[nr+2*i-1] + k2[i]*x[nr+2*i] + uu[nr+i];
            end for;

         else
            assert(false, "filterType (= " + String(filterType) + ") is unknown");
            uu = zeros(na+nr+1);
         end if;

         y = (gain*u_nominal)*uu[nr+na+1];

        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=true,
              extent={{-100.0,-100.0},{100.0,100.0}}),
              graphics={
            Line(points={{-80.0,80.0},{-80.0,-88.0}},
              color={192,192,192}),
            Polygon(lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid,
              points={{-80.0,92.0},{-88.0,70.0},{-72.0,70.0},{-80.0,92.0}}),
            Line(points={{-90.0,-78.0},{82.0,-78.0}},
              color={192,192,192}),
            Polygon(lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid,
              points={{90.0,-78.0},{68.0,-70.0},{68.0,-86.0},{90.0,-78.0}}),
            Text(textColor={192,192,192},
              extent={{-66.0,52.0},{88.0,90.0}},
              textString="%order"),
            Text(
              extent={{-138.0,-140.0},{162.0,-110.0}},
              textString="f_cut=%f_cut"),
            Rectangle(lineColor={160,160,164},
              fillColor={255,255,255},
              fillPattern=FillPattern.Backward,
              extent={{-80.0,-78.0},{22.0,10.0}}),
            Line(origin = {3.333,-6.667}, points = {{-83.333,34.667},{24.667,34.667},{42.667,-71.333}}, color = {0,0,127}, smooth = Smooth.Bezier)}),
          Documentation(info="<html>

<p>
This blocks models various types of filters:
</p>

<blockquote>
<strong>low pass, high pass, band pass, and band stop filters</strong>
</blockquote>

<p>
using various filter characteristics:
</p>

<blockquote>
<strong>CriticalDamping, Bessel, Butterworth, Chebyshev Type I filters</strong>
</blockquote>

<p>
By default, a filter block is initialized in <strong>steady-state</strong>, in order to
avoid unwanted oscillations at the beginning. In special cases, it might be
useful to select one of the other initialization options under tab
\"Advanced\".
</p>

<p>
Typical frequency responses for the 4 supported low pass filter types
are shown in the next figure:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Continuous/LowPassOrder4Filters.png\"
     alt=\"LowPassOrder4Filters.png\">
</blockquote>

<p>
The step responses of the same low pass filters are shown in the next figure,
starting from a steady state initial filter with initial input = 0.2:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Continuous/LowPassOrder4FiltersStepResponse.png\"
     alt=\"LowPassOrder4FiltersStepResponse.png\">
</blockquote>

<p>
Obviously, the frequency responses give a somewhat wrong impression
of the filter characteristics: Although Butterworth and Chebyshev
filters have a significantly steeper magnitude as the
CriticalDamping and Bessel filters, the step responses of
the latter ones are much better. This means for example, that
a CriticalDamping or a Bessel filter should be selected,
if a filter is mainly used to make a non-linear inverse model
realizable.
</p>

<p>
Typical frequency responses for the 4 supported high pass filter types
are shown in the next figure:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Continuous/HighPassOrder4Filters.png\"
     alt=\"HighPassOrder4Filters.png\">
</blockquote>

<p>
The corresponding step responses of these high pass filters are
shown in the next figure:
</p>
<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Continuous/HighPassOrder4FiltersStepResponse.png\"
     alt=\"HighPassOrder4FiltersStepResponse.png\">
</blockquote>

<p>
All filters are available in <strong>normalized</strong> (default) and non-normalized form.
In the normalized form, the amplitude of the filter transfer function
at the cut-off frequency f_cut is -3 dB (= 10^(-3/20) = 0.70794..).
Note, when comparing the filters of this function with other software systems,
the setting of \"normalized\" has to be selected appropriately. For example, the signal processing
toolbox of MATLAB provides the filters in non-normalized form and
therefore a comparison makes only sense, if normalized = <strong>false</strong>
is set. A normalized filter is usually better suited for applications,
since filters of different orders are \"comparable\",
whereas non-normalized filters usually require to adapt the
cut-off frequency, when the order of the filter is changed.
See a comparison of \"normalized\" and \"non-normalized\" filters at hand of
CriticalDamping filters of order 1,2,3:
</p>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Continuous/CriticalDampingNormalized.png\"
     alt=\"CriticalDampingNormalized.png\">
</blockquote>

<blockquote>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Continuous/CriticalDampingNonNormalized.png\"
     alt=\"CriticalDampingNonNormalized.png\">
</blockquote>

<h4>Implementation</h4>

<p>
The filters are implemented in the following, reliable way:
</p>

<ol>
<li> A prototype low pass filter with a cut-off angular frequency of 1 rad/s is constructed
     from the desired analogFilter and the desired normalization.</li>

<li> This prototype low pass filter is transformed to the desired filterType and the
     desired cut-off frequency f_cut using a transformation on the Laplace variable \"s\".</li>

<li> The resulting first and second order transfer functions are implemented in
     state space form, using the \"eigen value\" representation of a transfer function:
     <blockquote><pre>
// second order block with eigen values: a +/- jb
<strong>der</strong>(x1) = a*x1 - b*x2 + (a^2 + b^2)/b*u;
<strong>der</strong>(x2) = b*x1 + a*x2;
     y  = x2;
     </pre></blockquote>
     The dc-gain from the input to the output of this block is one and the selected
     states are in the order of the input (if \"u\" is in the order of \"one\", then the
     states are also in the order of \"one\"). In the \"Advanced\" tab, a \"nominal\" value for
     the input \"u\" can be given. If appropriately selected, the states are in the order of \"one\" and
     then step-size control is always appropriate.</li>
</ol>

<h4>References</h4>

<dl>
<dt>Tietze U., and Schenk C. (2002):</dt>
<dd> <strong>Halbleiter-Schaltungstechnik</strong>.
     Springer Verlag, 12. Auflage, pp. 815-852.</dd>
</dl>

</html>",     revisions="<html>
<dl>
  <dt><strong>Main Author:</strong></dt>
  <dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>,
      DLR Oberpfaffenhofen.</dd>
</dl>

<h4>Acknowledgement</h4>

<p>
The development of this block was partially funded by BMBF within the
     <a href=\"http://www.eurosyslib.com/\">ITEA2 EUROSYSLIB</a>
      project.
</p>

</html>"));
      end Filter;

      package Internal
      "Internal utility functions and blocks that should not be directly utilized by the user"
          extends Modelica.Icons.InternalPackage;

        package Filter
        "Internal utility functions for filters that should not be directly used"
            extends Modelica.Icons.InternalPackage;

          package base
          "Prototype low pass filters with cut-off frequency of 1 rad/s (other filters are derived by transformation from these base filters)"
              extends Modelica.Icons.InternalPackage;

          function CriticalDamping
              "Return base filter coefficients of CriticalDamping filter (= low pass filter with w_cut = 1 rad/s)"
            extends Modelica.Icons.Function;

            input Integer order(min=1) "Order of filter";
            input Boolean normalized=true
                "= true, if amplitude at f_cut = -3db, otherwise unmodified filter";

            output Real cr[order] "Coefficients of real poles";
            protected
            Real alpha=1.0 "Frequency correction factor";
            Real alpha2 "= alpha*alpha";
            Real den1[order]
                "[p] coefficients of denominator first order polynomials (a*p + 1)";
            Real den2[0,2]
                "[p^2, p] coefficients of denominator second order polynomials (b*p^2 + a*p + 1)";
            Real c0[0] "Coefficients of s^0 term if conjugate complex pole";
            Real c1[0] "Coefficients of s^1 term if conjugate complex pole";
          algorithm
            if normalized then
               // alpha := sqrt(2^(1/order) - 1);
               alpha := sqrt(10^(3/10/order)-1);
            else
               alpha := 1.0;
            end if;

            for i in 1:order loop
               den1[i] := alpha;
            end for;

            // Determine polynomials with highest power of s equal to one
              (cr,c0,c1) :=
                Modelica.Blocks.Continuous.Internal.Filter.Utilities.toHighestPowerOne(
                den1, den2);
          end CriticalDamping;

          function Bessel
              "Return base filter coefficients of Bessel filter (= low pass filter with w_cut = 1 rad/s)"
            extends Modelica.Icons.Function;

            input Integer order(min=1) "Order of filter";
            input Boolean normalized=true
                "= true, if amplitude at f_cut = -3db, otherwise unmodified filter";

            output Real cr[mod(order, 2)] "Coefficient of real pole";
            output Real c0[integer(order/2)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[integer(order/2)]
                "Coefficients of s^1 term if conjugate complex pole";
            protected
            Real alpha=1.0 "Frequency correction factor";
            Real alpha2 "= alpha*alpha";
            Real den1[size(cr,1)]
                "[p] coefficients of denominator first order polynomials (a*p + 1)";
            Real den2[size(c0, 1),2]
                "[p^2, p] coefficients of denominator second order polynomials (b*p^2 + a*p + 1)";
          algorithm
              (den1,den2,alpha) :=
                Modelica.Blocks.Continuous.Internal.Filter.Utilities.BesselBaseCoefficients(
                order);
            if not normalized then
               alpha2 := alpha*alpha;
               for i in 1:size(c0, 1) loop
                 den2[i, 1] := den2[i, 1]*alpha2;
                 den2[i, 2] := den2[i, 2]*alpha;
               end for;
               if size(cr,1) == 1 then
                 den1[1] := den1[1]*alpha;
               end if;
               end if;

            // Determine polynomials with highest power of s equal to one
              (cr,c0,c1) :=
                Modelica.Blocks.Continuous.Internal.Filter.Utilities.toHighestPowerOne(
                den1, den2);
          end Bessel;

          function Butterworth
              "Return base filter coefficients of Butterworth filter (= low pass filter with w_cut = 1 rad/s)"
            import Modelica.Constants.pi;
            extends Modelica.Icons.Function;

            input Integer order(min=1) "Order of filter";
            input Boolean normalized=true
                "= true, if amplitude at f_cut = -3db, otherwise unmodified filter";

            output Real cr[mod(order, 2)] "Coefficient of real pole";
            output Real c0[integer(order/2)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[integer(order/2)]
                "Coefficients of s^1 term if conjugate complex pole";
            protected
            Real alpha=1.0 "Frequency correction factor";
            Real alpha2 "= alpha*alpha";
            Real den1[size(cr,1)]
                "[p] coefficients of denominator first order polynomials (a*p + 1)";
            Real den2[size(c0, 1),2]
                "[p^2, p] coefficients of denominator second order polynomials (b*p^2 + a*p + 1)";
          algorithm
            for i in 1:size(c0, 1) loop
              den2[i, 1] := 1.0;
              den2[i, 2] := -2*Modelica.Math.cos(pi*(0.5 + (i - 0.5)/order));
            end for;
            if size(cr,1) == 1 then
              den1[1] := 1.0;
            end if;

            /* Transformation of filter transfer function with "new(p) = alpha*p"
     in order that the filter transfer function has an amplitude of
     -3 db at the cutoff frequency
  */
            /*
    if normalized then
      alpha := Internal.normalizationFactor(den1, den2);
      alpha2 := alpha*alpha;
      for i in 1:size(c0, 1) loop
        den2[i, 1] := den2[i, 1]*alpha2;
        den2[i, 2] := den2[i, 2]*alpha;
      end for;
      if size(cr,1) == 1 then
        den1[1] := den1[1]*alpha;
      end if;
    end if;
  */

            // Determine polynomials with highest power of s equal to one
              (cr,c0,c1) :=
                Modelica.Blocks.Continuous.Internal.Filter.Utilities.toHighestPowerOne(
                den1, den2);
          end Butterworth;

          function ChebyshevI
              "Return base filter coefficients of Chebyshev I filter (= low pass filter with w_cut = 1 rad/s)"
            import Modelica.Math.asinh;
            import Modelica.Constants.pi;

            extends Modelica.Icons.Function;

            input Integer order(min=1) "Order of filter";
            input Real A_ripple = 0.5 "Pass band ripple in [dB]";
            input Boolean normalized=true
                "= true, if amplitude at f_cut = -3db, otherwise unmodified filter";

            output Real cr[mod(order, 2)] "Coefficient of real pole";
            output Real c0[integer(order/2)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[integer(order/2)]
                "Coefficients of s^1 term if conjugate complex pole";
            protected
            Real epsilon;
            Real fac;
            Real alpha=1.0 "Frequency correction factor";
            Real alpha2 "= alpha*alpha";
            Real den1[size(cr,1)]
                "[p] coefficients of denominator first order polynomials (a*p + 1)";
            Real den2[size(c0, 1),2]
                "[p^2, p] coefficients of denominator second order polynomials (b*p^2 + a*p + 1)";
          algorithm
              epsilon := sqrt(10^(A_ripple/10) - 1);
              fac := asinh(1/epsilon)/order;

              den1 := fill(1/sinh(fac),size(den1,1));
              if size(cr,1) == 0 then
                 for i in 1:size(c0, 1) loop
                    den2[i,1] :=1/(cosh(fac)^2 - cos((2*i - 1)*pi/(2*order))^2);
                    den2[i,2] :=2*den2[i, 1]*sinh(fac)*cos((2*i - 1)*pi/(2*order));
                 end for;
              else
                 for i in 1:size(c0, 1) loop
                    den2[i,1] :=1/(cosh(fac)^2 - cos(i*pi/order)^2);
                    den2[i,2] :=2*den2[i, 1]*sinh(fac)*cos(i*pi/order);
                 end for;
              end if;

              /* Transformation of filter transfer function with "new(p) = alpha*p"
       in order that the filter transfer function has an amplitude of
       -3 db at the cutoff frequency
    */
              if normalized then
                alpha :=
                  Modelica.Blocks.Continuous.Internal.Filter.Utilities.normalizationFactor(
                  den1, den2);
                alpha2 := alpha*alpha;
                for i in 1:size(c0, 1) loop
                  den2[i, 1] := den2[i, 1]*alpha2;
                  den2[i, 2] := den2[i, 2]*alpha;
                end for;
                den1 := den1*alpha;
              end if;

            // Determine polynomials with highest power of s equal to one
              (cr,c0,c1) :=
                Modelica.Blocks.Continuous.Internal.Filter.Utilities.toHighestPowerOne(
                den1, den2);
          end ChebyshevI;
          end base;

          package coefficients "Filter coefficients"
              extends Modelica.Icons.InternalPackage;

          function lowPass
              "Return low pass filter coefficients at given cut-off frequency"
            import Modelica.Constants.pi;
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles";
            input Real c0_in[:]
                "Coefficients of s^0 term if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";
            input SI.Frequency f_cut "Cut-off frequency";

            output Real cr[size(cr_in,1)] "Coefficient of real pole";
            output Real c0[size(c0_in,1)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";

            protected
            SI.AngularVelocity w_cut=2*pi*f_cut
                "Cut-off angular frequency";
            Real w_cut2=w_cut*w_cut;

          algorithm
            assert(f_cut > 0, "Cut-off frequency f_cut must be positive");

            /* Change filter coefficients according to transformation new(s) = s/w_cut
     s + cr           -> (s/w) + cr              = (s + w*cr)/w
     s^2 + c1*s + c0  -> (s/w)^2 + c1*(s/w) + c0 = (s^2 + (c1*w)*s + (c0*w^2))/w^2
  */
            cr := w_cut*cr_in;
            c1 := w_cut*c1_in;
            c0 := w_cut2*c0_in;

          end lowPass;

          function highPass
              "Return high pass filter coefficients at given cut-off frequency"
            import Modelica.Constants.pi;
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles";
            input Real c0_in[:]
                "Coefficients of s^0 term if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";
            input SI.Frequency f_cut "Cut-off frequency";

            output Real cr[size(cr_in,1)] "Coefficient of real pole";
            output Real c0[size(c0_in,1)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";

            protected
            SI.AngularVelocity w_cut=2*pi*f_cut
                "Cut-off angular frequency";
            Real w_cut2=w_cut*w_cut;

          algorithm
            assert(f_cut > 0, "Cut-off frequency f_cut must be positive");

            /* Change filter coefficients according to transformation: new(s) = 1/s
        1/(s + cr)          -> 1/(1/s + cr)                = (1/cr)*s / (s + (1/cr))
        1/(s^2 + c1*s + c0) -> 1/((1/s)^2 + c1*(1/s) + c0) = (1/c0)*s^2 / (s^2 + (c1/c0)*s + 1/c0)

     Check whether transformed roots are also conjugate complex:
        c0 - c1^2/4 > 0  -> (1/c0) - (c1/c0)^2 / 4
                            = (c0 - c1^2/4) / c0^2 > 0
        It is therefore guaranteed that the roots remain conjugate complex

     Change filter coefficients according to transformation new(s) = s/w_cut
        s + 1/cr                -> (s/w) + 1/cr                   = (s + w/cr)/w
        s^2 + (c1/c0)*s + 1/c0  -> (s/w)^2 + (c1/c0)*(s/w) + 1/c0 = (s^2 + (w*c1/c0)*s + (w^2/c0))/w^2
  */
            for i in 1:size(cr_in,1) loop
               cr[i] := w_cut/cr_in[i];
            end for;

            for i in 1:size(c0_in,1) loop
               c0[i] := w_cut2/c0_in[i];
               c1[i] := w_cut*c1_in[i]/c0_in[i];
            end for;

          end highPass;

          function bandPass
              "Return band pass filter coefficients at given cut-off frequency"
            import Modelica.Constants.pi;
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles";
            input Real c0_in[:]
                "Coefficients of s^0 term if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";
            input SI.Frequency f_min
                "Band of band pass filter is f_min (A=-3db) .. f_max (A=-3db)";
            input SI.Frequency f_max "Upper band frequency";

            output Real cr[0] "Coefficient of real pole";
            output Real c0[size(cr_in,1) + 2*size(c0_in,1)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[size(cr_in,1) + 2*size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";
            output Real cn "Numerator coefficient of the PT2 terms";
            protected
            SI.Frequency f0 = sqrt(f_min*f_max);
            SI.AngularVelocity w_cut=2*pi*f0
                "Cut-off angular frequency";
            Real w_band = (f_max - f_min) / f0;
            Real w_cut2=w_cut*w_cut;
            Real c;
            Real alpha;
            Integer j;
          algorithm
            assert(f_min > 0 and f_min < f_max, "Band frequencies f_min and f_max are wrong");

              /* The band pass filter is derived from the low pass filter by
       the transformation new(s) = (s + 1/s)/w   (w = w_band = (f_max - f_min)/sqrt(f_max*f_min) )

       1/(s + cr)         -> 1/((s/w + 1/s/w) + cr)
                             = w*s / (s^2 + cr*w*s + 1)

       1/(s^2 + c1*s + c0) -> 1/( (s+1/s)^2/w^2 + c1*(s + 1/s)/w + c0 )
                              = 1 /( ( s^2 + 1/s^2 + 2)/w^2 + (s + 1/s)*c1/w + c0 )
                              = w^2*s^2 / (s^4 + 2*s^2 + 1 + (s^3 + s)*c1*w + c0*w^2*s^2)
                              = w^2*s^2 / (s^4 + c1*w*s^3 + (2+c0*w^2)*s^2 + c1*w*s + 1)

                              Assume the following description with PT2:
                              = w^2*s^2 /( (s^2 + s*(c/alpha) + 1/alpha^2)*
                                           (s^2 + s*(c*alpha) + alpha^2) )
                              = w^2*s^2 / ( s^4 + c*(alpha + 1/alpha)*s^3
                                                + (alpha^2 + 1/alpha^2 + c^2)*s^2
                                                + c*(alpha + 1/alpha)*s + 1 )

                              and therefore:
                                c*(alpha + 1/alpha) = c1*w       -> c = c1*w / (alpha + 1/alpha)
                                                                      = c1*w*alpha/(1+alpha^2)
                                alpha^2 + 1/alpha^2 + c^2 = 2+c0*w^2 -> equation to determine alpha
                                alpha^4 + 1 + c1^2*w^2*alpha^4/(1+alpha^2)^2 = (2+c0*w^2)*alpha^2
                                or z = alpha^2
                                z^2 + c^1^2*w^2*z^2/(1+z)^2 - (2+c0*w^2)*z + 1 = 0

     Check whether roots remain conjugate complex
        c0 - (c1/2)^2 > 0:    1/alpha^2 - (c/alpha)^2/4
                              = 1/alpha^2*(1 - c^2/4)    -> not possible to figure this out

     Afterwards, change filter coefficients according to transformation new(s) = s/w_cut
        w_band*s/(s^2 + c1*s + c0)  -> w_band*(s/w)/((s/w)^2 + c1*(s/w) + c0 =
                                       (w_band/w)*s/(s^2 + (c1*w)*s + (c0*w^2))/w^2) =
                                       (w_band*w)*s/(s^2 + (c1*w)*s + (c0*w^2))
    */
              for i in 1:size(cr_in,1) loop
                 c1[i] := w_cut*cr_in[i]*w_band;
                 c0[i] := w_cut2;
              end for;

              for i in 1:size(c1_in,1) loop
                alpha :=
                  Modelica.Blocks.Continuous.Internal.Filter.Utilities.bandPassAlpha(
                        c1_in[i],
                        c0_in[i],
                        w_band);
                 c       := c1_in[i]*w_band / (alpha + 1/alpha);
                 j       := size(cr_in,1) + 2*i - 1;
                 c1[j]   := w_cut*c/alpha;
                 c1[j+1] := w_cut*c*alpha;
                 c0[j]   := w_cut2/alpha^2;
                 c0[j+1] := w_cut2*alpha^2;
              end for;

              cn :=w_band*w_cut;

          end bandPass;

          function bandStop
              "Return band stop filter coefficients at given cut-off frequency"
            import Modelica.Constants.pi;
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles";
            input Real c0_in[:]
                "Coefficients of s^0 term if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";
            input SI.Frequency f_min
                "Band of band stop filter is f_min (A=-3db) .. f_max (A=-3db)";
            input SI.Frequency f_max "Upper band frequency";

            output Real cr[0] "Coefficient of real pole";
            output Real c0[size(cr_in,1) + 2*size(c0_in,1)]
                "Coefficients of s^0 term if conjugate complex pole";
            output Real c1[size(cr_in,1) + 2*size(c0_in,1)]
                "Coefficients of s^1 term if conjugate complex pole";
            protected
            SI.Frequency f0 = sqrt(f_min*f_max);
            SI.AngularVelocity w_cut=2*pi*f0
                "Cut-off angular frequency";
            Real w_band = (f_max - f_min) / f0;
            Real w_cut2=w_cut*w_cut;
            Real c;
            Real ww;
            Real alpha;
            Integer j;
          algorithm
            assert(f_min > 0 and f_min < f_max, "Band frequencies f_min and f_max are wrong");

              /* The band pass filter is derived from the low pass filter by
       the transformation new(s) = (s + 1/s)/w   (w = w_band = (f_max - f_min)/sqrt(f_max*f_min) )

       1/(s + cr)         -> 1/((s/w + 1/s/w) + cr)
                             = w*s / (s^2 + cr*w*s + 1)

       1/(s^2 + c1*s + c0) -> 1/( (s+1/s)^2/w^2 + c1*(s + 1/s)/w + c0 )
                              = 1 /( ( s^2 + 1/s^2 + 2)/w^2 + (s + 1/s)*c1/w + c0 )
                              = w^2*s^2 / (s^4 + 2*s^2 + 1 + (s^3 + s)*c1*w + c0*w^2*s^2)
                              = w^2*s^2 / (s^4 + c1*w*s^3 + (2+c0*w^2)*s^2 + c1*w*s + 1)

                              Assume the following description with PT2:
                              = w^2*s^2 /( (s^2 + s*(c/alpha) + 1/alpha^2)*
                                           (s^2 + s*(c*alpha) + alpha^2) )
                              = w^2*s^2 / ( s^4 + c*(alpha + 1/alpha)*s^3
                                                + (alpha^2 + 1/alpha^2 + c^2)*s^2
                                                + c*(alpha + 1/alpha)*s + 1 )

                              and therefore:
                                c*(alpha + 1/alpha) = c1*w       -> c = c1*w / (alpha + 1/alpha)
                                                                      = c1*w*alpha/(1+alpha^2)
                                alpha^2 + 1/alpha^2 + c^2 = 2+c0*w^2 -> equation to determine alpha
                                alpha^4 + 1 + c1^2*w^2*alpha^4/(1+alpha^2)^2 = (2+c0*w^2)*alpha^2
                                or z = alpha^2
                                z^2 + c^1^2*w^2*z^2/(1+z)^2 - (2+c0*w^2)*z + 1 = 0

       The band stop filter is derived from the low pass filter by
       the transformation new(s) = w/( (s + 1/s) )   (w = w_band = (f_max - f_min)/sqrt(f_max*f_min) )

       cr/(s + cr)         -> 1/(( w/(s + 1/s) ) + cr)
                              = (s^2 + 1) / (s^2 + (w/cr)*s + 1)

       c0/(s^2 + c1*s + c0) -> c0/( w^2/(s + 1/s)^2 + c1*w/(s + 1/s) + c0 )
                               = c0*(s^2 + 1)^2 / (s^4 + c1*w*s^3/c0 + (2+w^2/b)*s^2 + c1*w*s/c0 + 1)

                               Assume the following description with PT2:
                               = c0*(s^2 + 1)^2 / ( (s^2 + s*(c/alpha) + 1/alpha^2)*
                                                    (s^2 + s*(c*alpha) + alpha^2) )
                               = c0*(s^2 + 1)^2 / (  s^4 + c*(alpha + 1/alpha)*s^3
                                                         + (alpha^2 + 1/alpha^2 + c^2)*s^2
                                                         + c*(alpha + 1/alpha)*p + 1 )

                            and therefore:
                              c*(alpha + 1/alpha) = c1*w/b         -> c = c1*w/(c0*(alpha + 1/alpha))
                              alpha^2 + 1/alpha^2 + c^2 = 2+w^2/c0 -> equation to determine alpha
                              alpha^4 + 1 + (c1*w/c0*alpha^2)^2/(1+alpha^2)^2 = (2+w^2/c0)*alpha^2
                              or z = alpha^2
                              z^2 + (c1*w/c0*z)^2/(1+z)^2 - (2+w^2/c0)*z + 1 = 0

                            same as:  ww = w/c0
                              z^2 + (c1*ww*z)^2/(1+z)^2 - (2+c0*ww)*z + 1 = 0  -> same equation as for BandPass

     Afterwards, change filter coefficients according to transformation new(s) = s/w_cut
        c0*(s^2+1)(s^2 + c1*s + c0)  -> c0*((s/w)^2 + 1) / ((s/w)^2 + c1*(s/w) + c0 =
                                        c0/w^2*(s^2 + w^2) / (s^2 + (c1*w)*s + (c0*w^2))/w^2) =
                                        (s^2 + c0*w^2) / (s^2 + (c1*w)*s + (c0*w^2))
    */
              for i in 1:size(cr_in,1) loop
                 c1[i] := w_cut*w_band/cr_in[i];
                 c0[i] := w_cut2;
              end for;

              for i in 1:size(c1_in,1) loop
                 ww      := w_band/c0_in[i];
                alpha :=
                  Modelica.Blocks.Continuous.Internal.Filter.Utilities.bandPassAlpha(
                        c1_in[i],
                        c0_in[i],
                        ww);
                 c       := c1_in[i]*ww / (alpha + 1/alpha);
                 j       := size(cr_in,1) + 2*i - 1;
                 c1[j]   := w_cut*c/alpha;
                 c1[j+1] := w_cut*c*alpha;
                 c0[j]   := w_cut2/alpha^2;
                 c0[j+1] := w_cut2*alpha^2;
              end for;

          end bandStop;
          end coefficients;

          package roots
          "Filter roots and gain as needed for block implementations"
              extends Modelica.Icons.InternalPackage;

          function lowPass
              "Return low pass filter roots as needed for block for given cut-off frequency"
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles of base filter";
            input Real c0_in[:]
                "Coefficients of s^0 term of base filter if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term of base filter if conjugate complex pole";
            input SI.Frequency f_cut "Cut-off frequency";

            output Real r[size(cr_in,1)] "Real eigenvalues";
            output Real a[size(c0_in,1)]
                "Real parts of complex conjugate eigenvalues";
            output Real b[size(c0_in,1)]
                "Imaginary parts of complex conjugate eigenvalues";
            output Real ku[size(c0_in,1)] "Input gain";
            protected
            Real c0[size(c0_in,1)];
            Real c1[size(c0_in,1)];
            Real cr[size(cr_in,1)];
          algorithm
            // Get coefficients of low pass filter at f_cut
            (cr, c0, c1) :=coefficients.lowPass(cr_in, c0_in, c1_in, f_cut);

            // Transform coefficients in to root
            for i in 1:size(cr_in,1) loop
              r[i] :=-cr[i];
            end for;

            for i in 1:size(c0_in,1) loop
              a [i] :=-c1[i]/2;
              b [i] :=sqrt(c0[i] - a[i]*a[i]);
              ku[i] :=c0[i]/b[i];
            end for;

            annotation (Documentation(info="<html>

<p>
The goal is to implement the filter in the following form:
</p>

<blockquote><pre>
// real pole:
 der(x) = r*x - r*u
     y  = x

// complex conjugate poles:
der(x1) = a*x1 - b*x2 + ku*u;
der(x2) = b*x1 + a*x2;
     y  = x2;

          ku = (a^2 + b^2)/b
</pre></blockquote>
<p>
This representation has the following transfer function:
</p>
<blockquote><pre>
// real pole:
    s*y = r*y - r*u
  or
    (s-r)*y = -r*u
  or
    y = -r/(s-r)*u

  comparing coefficients with
    y = cr/(s + cr)*u  ->  r = -cr      // r is the real eigenvalue

// complex conjugate poles
    s*x2 =  a*x2 + b*x1
    s*x1 = -b*x2 + a*x1 + ku*u
  or
    (s-a)*x2               = b*x1  ->  x2 = b/(s-a)*x1
    (s + b^2/(s-a) - a)*x1 = ku*u  ->  (s(s-a) + b^2 - a*(s-a))*x1  = ku*(s-a)*u
                                   ->  (s^2 - 2*a*s + a^2 + b^2)*x1 = ku*(s-a)*u
  or
    x1 = ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
    x2 = b/(s-a)*ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
       = b*ku/(s^2 - 2*a*s + a^2 + b^2)*u
    y  = x2

  comparing coefficients with
    y = c0/(s^2 + c1*s + c0)*u  ->  a  = -c1/2
                                    b  = sqrt(c0 - a^2)
                                    ku = c0/b
                                       = (a^2 + b^2)/b

  comparing with eigenvalue representation:
    (s - (a+jb))*(s - (a-jb)) = s^2 -2*a*s + a^2 + b^2
  shows that:
    a: real part of eigenvalue
    b: imaginary part of eigenvalue

  time -> infinity:
    y(s=0) = x2(s=0) = 1
             x1(s=0) = -ku*a/(a^2 + b^2)*u
                     = -(a/b)*u
</pre></blockquote>

</html>"));
          end lowPass;

          function highPass
              "Return high pass filter roots as needed for block for given cut-off frequency"
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles of base filter";
            input Real c0_in[:]
                "Coefficients of s^0 term of base filter if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term of base filter if conjugate complex pole";
            input SI.Frequency f_cut "Cut-off frequency";

            output Real r[size(cr_in,1)] "Real eigenvalues";
            output Real a[size(c0_in,1)]
                "Real parts of complex conjugate eigenvalues";
            output Real b[size(c0_in,1)]
                "Imaginary parts of complex conjugate eigenvalues";
            output Real ku[size(c0_in,1)] "Gains of input terms";
            output Real k1[size(c0_in,1)] "Gains of y = k1*x1 + k2*x + u";
            output Real k2[size(c0_in,1)] "Gains of y = k1*x1 + k2*x + u";
            protected
            Real c0[size(c0_in,1)];
            Real c1[size(c0_in,1)];
            Real cr[size(cr_in,1)];
            Real ba2;
          algorithm
            // Get coefficients of high pass filter at f_cut
            (cr, c0, c1) :=coefficients.highPass(cr_in, c0_in, c1_in, f_cut);

            // Transform coefficients in to roots
            for i in 1:size(cr_in,1) loop
              r[i] :=-cr[i];
            end for;

            for i in 1:size(c0_in,1) loop
              a[i]  := -c1[i]/2;
              b[i]  := sqrt(c0[i] - a[i]*a[i]);
              ku[i] := c0[i]/b[i];
              k1[i] := 2*a[i]/ku[i];
              ba2   := (b[i]/a[i])^2;
              k2[i] := (1-ba2)/(1+ba2);
            end for;

            annotation (Documentation(info="<html>

<p>
The goal is to implement the filter in the following form:
</p>

<blockquote><pre>
// real pole:
 der(x) = r*x - r*u
     y  = -x + u

// complex conjugate poles:
der(x1) = a*x1 - b*x2 + ku*u;
der(x2) = b*x1 + a*x2;
     y  = k1*x1 + k2*x2 + u;

          ku = (a^2 + b^2)/b
          k1 = 2*a/ku
          k2 = (a^2 - b^2) / (b*ku)
             = (a^2 - b^2) / (a^2 + b^2)
             = (1 - (b/a)^2) / (1 + (b/a)^2)
</pre></blockquote>
<p>
This representation has the following transfer function:
</p>
<blockquote><pre>
// real pole:
    s*x = r*x - r*u
  or
    (s-r)*x = -r*u   -> x = -r/(s-r)*u
  or
    y = r/(s-r)*u + (s-r)/(s-r)*u
      = (r+s-r)/(s-r)*u
      = s/(s-r)*u

// comparing coefficients with
    y = s/(s + cr)*u  ->  r = -cr      // r is the real eigenvalue

// complex conjugate poles
    s*x2 =  a*x2 + b*x1
    s*x1 = -b*x2 + a*x1 + ku*u
  or
    (s-a)*x2               = b*x1  ->  x2 = b/(s-a)*x1
    (s + b^2/(s-a) - a)*x1 = ku*u  ->  (s(s-a) + b^2 - a*(s-a))*x1  = ku*(s-a)*u
                                   ->  (s^2 - 2*a*s + a^2 + b^2)*x1 = ku*(s-a)*u
  or
    x1 = ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
    x2 = b/(s-a)*ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
       = b*ku/(s^2 - 2*a*s + a^2 + b^2)*u
    y  = k1*x1 + k2*x2 + u
       = (k1*ku*(s-a) + k2*b*ku +  s^2 - 2*a*s + a^2 + b^2) /
         (s^2 - 2*a*s + a^2 + b^2)*u
       = (s^2 + (k1*ku - 2*a)*s + k2*b*ku - k1*ku*a + a^2 + b^2) /
         (s^2 - 2*a*s + a^2 + b^2)*u
       = (s^2 + (2*a-2*a)*s + a^2 - b^2 - 2*a^2 + a^2 + b^2) /
         (s^2 - 2*a*s + a^2 + b^2)*u
       = s^2 / (s^2 - 2*a*s + a^2 + b^2)*u

// comparing coefficients with
    y = s^2/(s^2 + c1*s + c0)*u  ->  a = -c1/2
                                     b = sqrt(c0 - a^2)

// comparing with eigenvalue representation:
    (s - (a+jb))*(s - (a-jb)) = s^2 -2*a*s + a^2 + b^2
// shows that:
//   a: real part of eigenvalue
//   b: imaginary part of eigenvalue
</pre></blockquote>

</html>"));
          end highPass;

          function bandPass
              "Return band pass filter roots as needed for block for given cut-off frequency"
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles of base filter";
            input Real c0_in[:]
                "Coefficients of s^0 term of base filter if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term of base filter if conjugate complex pole";
            input SI.Frequency f_min
                "Band of band pass filter is f_min (A=-3db) .. f_max (A=-3db)";
            input SI.Frequency f_max "Upper band frequency";

            output Real a[size(cr_in,1) + 2*size(c0_in,1)]
                "Real parts of complex conjugate eigenvalues";
            output Real b[size(cr_in,1) + 2*size(c0_in,1)]
                "Imaginary parts of complex conjugate eigenvalues";
            output Real ku[size(cr_in,1) + 2*size(c0_in,1)] "Gains of input terms";
            output Real k1[size(cr_in,1) + 2*size(c0_in,1)]
                "Gains of y = k1*x1 + k2*x";
            output Real k2[size(cr_in,1) + 2*size(c0_in,1)]
                "Gains of y = k1*x1 + k2*x";
            protected
            Real cr[0];
            Real c0[size(a,1)];
            Real c1[size(a,1)];
            Real cn;
            Real bb;
          algorithm
            // Get coefficients of band pass filter at f_cut
            (cr, c0, c1, cn) :=coefficients.bandPass(cr_in, c0_in, c1_in, f_min, f_max);

            // Transform coefficients in to roots
            for i in 1:size(a,1) loop
              a[i]  := -c1[i]/2;
              bb    := c0[i] - a[i]*a[i];
              assert(bb >= 0, "\nNot possible to use band pass filter, since transformation results in\n"+
                              "system that does not have conjugate complex poles.\n" +
                              "Try to use another analog filter for the band pass.\n");
              b[i]  := sqrt(bb);
              ku[i] := c0[i]/b[i];
              k1[i] := cn/ku[i];
              k2[i] := cn*a[i]/(b[i]*ku[i]);
            end for;

            annotation (Documentation(info="<html>

<p>
The goal is to implement the filter in the following form:
</p>

<blockquote><pre>
// complex conjugate poles:
der(x1) = a*x1 - b*x2 + ku*u;
der(x2) = b*x1 + a*x2;
     y  = k1*x1 + k2*x2;

          ku = (a^2 + b^2)/b
          k1 = cn/ku
          k2 = cn*a/(b*ku)
</pre></blockquote>
<p>
This representation has the following transfer function:
</p>
<blockquote><pre>
// complex conjugate poles
    s*x2 =  a*x2 + b*x1
    s*x1 = -b*x2 + a*x1 + ku*u
  or
    (s-a)*x2               = b*x1  ->  x2 = b/(s-a)*x1
    (s + b^2/(s-a) - a)*x1 = ku*u  ->  (s(s-a) + b^2 - a*(s-a))*x1  = ku*(s-a)*u
                                   ->  (s^2 - 2*a*s + a^2 + b^2)*x1 = ku*(s-a)*u
  or
    x1 = ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
    x2 = b/(s-a)*ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
       = b*ku/(s^2 - 2*a*s + a^2 + b^2)*u
    y  = k1*x1 + k2*x2
       = (k1*ku*(s-a) + k2*b*ku) / (s^2 - 2*a*s + a^2 + b^2)*u
       = (k1*ku*s + k2*b*ku - k1*ku*a) / (s^2 - 2*a*s + a^2 + b^2)*u
       = (cn*s + cn*a - cn*a) / (s^2 - 2*a*s + a^2 + b^2)*u
       = cn*s / (s^2 - 2*a*s + a^2 + b^2)*u

  comparing coefficients with
    y = cn*s / (s^2 + c1*s + c0)*u  ->  a = -c1/2
                                        b = sqrt(c0 - a^2)

  comparing with eigenvalue representation:
    (s - (a+jb))*(s - (a-jb)) = s^2 -2*a*s + a^2 + b^2
  shows that:
    a: real part of eigenvalue
    b: imaginary part of eigenvalue
</pre></blockquote>

</html>"));
          end bandPass;

          function bandStop
              "Return band stop filter roots as needed for block for given cut-off frequency"
            extends Modelica.Icons.Function;

            input Real cr_in[:] "Coefficients of real poles of base filter";
            input Real c0_in[:]
                "Coefficients of s^0 term of base filter if conjugate complex pole";
            input Real c1_in[size(c0_in,1)]
                "Coefficients of s^1 term of base filter if conjugate complex pole";
            input SI.Frequency f_min
                "Band of band stop filter is f_min (A=-3db) .. f_max (A=-3db)";
            input SI.Frequency f_max "Upper band frequency";

            output Real a[size(cr_in,1) + 2*size(c0_in,1)]
                "Real parts of complex conjugate eigenvalues";
            output Real b[size(cr_in,1) + 2*size(c0_in,1)]
                "Imaginary parts of complex conjugate eigenvalues";
            output Real ku[size(cr_in,1) + 2*size(c0_in,1)] "Gains of input terms";
            output Real k1[size(cr_in,1) + 2*size(c0_in,1)]
                "Gains of y = k1*x1 + k2*x";
            output Real k2[size(cr_in,1) + 2*size(c0_in,1)]
                "Gains of y = k1*x1 + k2*x";
            protected
            Real cr[0];
            Real c0[size(a,1)];
            Real c1[size(a,1)];
            Real cn;
            Real bb;
          algorithm
            // Get coefficients of band stop filter at f_cut
            (cr, c0, c1) :=coefficients.bandStop(cr_in, c0_in, c1_in, f_min, f_max);

            // Transform coefficients in to roots
            for i in 1:size(a,1) loop
              a[i]  := -c1[i]/2;
              bb    := c0[i] - a[i]*a[i];
              assert(bb >= 0, "\nNot possible to use band stop filter, since transformation results in\n"+
                              "system that does not have conjugate complex poles.\n" +
                              "Try to use another analog filter for the band stop filter.\n");
              b[i]  := sqrt(bb);
              ku[i] := c0[i]/b[i];
              k1[i] := 2*a[i]/ku[i];
              k2[i] := (c0[i] + a[i]^2 - b[i]^2)/(b[i]*ku[i]);
            end for;

            annotation (Documentation(info="<html>

<p>
The goal is to implement the filter in the following form:
</p>

<blockquote><pre>
// complex conjugate poles:
der(x1) = a*x1 - b*x2 + ku*u;
der(x2) = b*x1 + a*x2;
     y  = k1*x1 + k2*x2 + u;

          ku = (a^2 + b^2)/b
          k1 = 2*a/ku
          k2 = (c0 + a^2 - b^2)/(b*ku)
</pre></blockquote>
<p>
This representation has the following transfer function:
</p>
<blockquote><pre>
// complex conjugate poles
    s*x2 =  a*x2 + b*x1
    s*x1 = -b*x2 + a*x1 + ku*u
  or
    (s-a)*x2               = b*x1  ->  x2 = b/(s-a)*x1
    (s + b^2/(s-a) - a)*x1 = ku*u  ->  (s(s-a) + b^2 - a*(s-a))*x1  = ku*(s-a)*u
                                   ->  (s^2 - 2*a*s + a^2 + b^2)*x1 = ku*(s-a)*u
  or
    x1 = ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
    x2 = b/(s-a)*ku*(s-a)/(s^2 - 2*a*s + a^2 + b^2)*u
       = b*ku/(s^2 - 2*a*s + a^2 + b^2)*u
    y  = k1*x1 + k2*x2 + u
       = (k1*ku*(s-a) + k2*b*ku + s^2 - 2*a*s + a^2 + b^2) / (s^2 - 2*a*s + a^2 + b^2)*u
       = (s^2 + (k1*ku-2*a)*s + k2*b*ku - k1*ku*a + a^2 + b^2) / (s^2 - 2*a*s + a^2 + b^2)*u
       = (s^2 + c0 + a^2 - b^2 - 2*a^2 + a^2 + b^2) / (s^2 - 2*a*s + a^2 + b^2)*u
       = (s^2 + c0) / (s^2 - 2*a*s + a^2 + b^2)*u

  comparing coefficients with
    y = (s^2 + c0) / (s^2 + c1*s + c0)*u  ->  a = -c1/2
                                              b = sqrt(c0 - a^2)

  comparing with eigenvalue representation:
    (s - (a+jb))*(s - (a-jb)) = s^2 -2*a*s + a^2 + b^2
  shows that:
    a: real part of eigenvalue
    b: imaginary part of eigenvalue
</pre></blockquote>

</html>"));
          end bandStop;
          end roots;

          package Utilities "Utility functions for filter computations"
              extends Modelica.Icons.InternalPackage;

            function BesselBaseCoefficients
              "Return coefficients of normalized low pass Bessel filter (= gain at cut-off frequency 1 rad/s is decreased 3dB)"
              extends Modelica.Icons.Function;

              import Modelica.Utilities.Streams;
              input Integer order "Order of filter in the range 1..41";
              output Real c1[mod(order, 2)]
                "[p] coefficients of Bessel denominator polynomials (a*p + 1)";
              output Real c2[integer(order/2),2]
                "[p^2, p] coefficients of Bessel denominator polynomials (b2*p^2 + b1*p + 1)";
              output Real alpha "Normalization factor";
            algorithm
              if order == 1 then
                alpha := 1.002377293007601;
                c1[1] := 0.9976283451109835;
              elseif order == 2 then
                alpha := 0.7356641785819585;
                c2[1, 1] := 0.6159132201783791;
                c2[1, 2] := 1.359315879600889;
              elseif order == 3 then
                alpha := 0.5704770156982642;
                c1[1] := 0.7548574865985343;
                c2[1, 1] := 0.4756958028827457;
                c2[1, 2] := 0.9980615136104388;
              elseif order == 4 then
                alpha := 0.4737978580281427;
                c2[1, 1] := 0.4873729247240677;
                c2[1, 2] := 1.337564170455762;
                c2[2, 1] := 0.3877724315741958;
                c2[2, 2] := 0.7730405590839861;
              elseif order == 5 then
                alpha := 0.4126226974763408;
                c1[1] := 0.6645723262620757;
                c2[1, 1] := 0.4115231900614016;
                c2[1, 2] := 1.138349926728708;
                c2[2, 1] := 0.3234938702877912;
                c2[2, 2] := 0.6205992985771313;
              elseif order == 6 then
                alpha := 0.3705098000736233;
                c2[1, 1] := 0.3874508649098960;
                c2[1, 2] := 1.219740879520741;
                c2[2, 1] := 0.3493298843155746;
                c2[2, 2] := 0.9670265529381365;
                c2[3, 1] := 0.2747419229514599;
                c2[3, 2] := 0.5122165075105700;
              elseif order == 7 then
                alpha := 0.3393452623586350;
                c1[1] := 0.5927147125821412;
                c2[1, 1] := 0.3383379423919174;
                c2[1, 2] := 1.092630816438030;
                c2[2, 1] := 0.3001025788696046;
                c2[2, 2] := 0.8289928256598656;
                c2[3, 1] := 0.2372867471539579;
                c2[3, 2] := 0.4325128641920154;
              elseif order == 8 then
                alpha := 0.3150267393795002;
                c2[1, 1] := 0.3151115975207653;
                c2[1, 2] := 1.109403015460190;
                c2[2, 1] := 0.2969344839572762;
                c2[2, 2] := 0.9737455812222699;
                c2[3, 1] := 0.2612545921889538;
                c2[3, 2] := 0.7190394712068573;
                c2[4, 1] := 0.2080523342974281;
                c2[4, 2] := 0.3721456473047434;
              elseif order == 9 then
                alpha := 0.2953310177184124;
                c1[1] := 0.5377196679501422;
                c2[1, 1] := 0.2824689124281034;
                c2[1, 2] := 1.022646191567475;
                c2[2, 1] := 0.2626824161383468;
                c2[2, 2] := 0.8695626454762596;
                c2[3, 1] := 0.2302781917677917;
                c2[3, 2] := 0.6309047553448520;
                c2[4, 1] := 0.1847991729757028;
                c2[4, 2] := 0.3251978031287202;
              elseif order == 10 then
                alpha := 0.2789426890619463;
                c2[1, 1] := 0.2640769908255582;
                c2[1, 2] := 1.019788132875305;
                c2[2, 1] := 0.2540802639216947;
                c2[2, 2] := 0.9377020417760623;
                c2[3, 1] := 0.2343577229427963;
                c2[3, 2] := 0.7802229808216112;
                c2[4, 1] := 0.2052193139338624;
                c2[4, 2] := 0.5594176813008133;
                c2[5, 1] := 0.1659546953748916;
                c2[5, 2] := 0.2878349616233292;
              elseif order == 11 then
                alpha := 0.2650227766037203;
                c1[1] := 0.4950265498954191;
                c2[1, 1] := 0.2411858478546218;
                c2[1, 2] := 0.9567800996387417;
                c2[2, 1] := 0.2296849355380925;
                c2[2, 2] := 0.8592523717113126;
                c2[3, 1] := 0.2107851705677406;
                c2[3, 2] := 0.7040216048898129;
                c2[4, 1] := 0.1846461385164021;
                c2[4, 2] := 0.5006729207276717;
                c2[5, 1] := 0.1504217970817433;
                c2[5, 2] := 0.2575070491320295;
              elseif order == 12 then
                alpha := 0.2530051198547209;
                c2[1, 1] := 0.2268294941204543;
                c2[1, 2] := 0.9473116570034053;
                c2[2, 1] := 0.2207657387793729;
                c2[2, 2] := 0.8933728946287606;
                c2[3, 1] := 0.2087600700376653;
                c2[3, 2] := 0.7886236252756229;
                c2[4, 1] := 0.1909959101492760;
                c2[4, 2] := 0.6389263649257017;
                c2[5, 1] := 0.1675208146048472;
                c2[5, 2] := 0.4517847275162215;
                c2[6, 1] := 0.1374257286372761;
                c2[6, 2] := 0.2324699157474680;
              elseif order == 13 then
                alpha := 0.2424910397561007;
                c1[1] := 0.4608848369928040;
                c2[1, 1] := 0.2099813050274780;
                c2[1, 2] := 0.8992478823790660;
                c2[2, 1] := 0.2027250423101359;
                c2[2, 2] := 0.8328117484224146;
                c2[3, 1] := 0.1907635894058731;
                c2[3, 2] := 0.7257379204691213;
                c2[4, 1] := 0.1742280397887686;
                c2[4, 2] := 0.5830640944868014;
                c2[5, 1] := 0.1530858190490478;
                c2[5, 2] := 0.4106192089751885;
                c2[6, 1] := 0.1264090712880446;
                c2[6, 2] := 0.2114980230156001;
              elseif order == 14 then
                alpha := 0.2331902368695848;
                c2[1, 1] := 0.1986162311411235;
                c2[1, 2] := 0.8876961808055535;
                c2[2, 1] := 0.1946683341271615;
                c2[2, 2] := 0.8500754229171967;
                c2[3, 1] := 0.1868331332895056;
                c2[3, 2] := 0.7764629313723603;
                c2[4, 1] := 0.1752118757862992;
                c2[4, 2] := 0.6699720402924552;
                c2[5, 1] := 0.1598906457908402;
                c2[5, 2] := 0.5348446712848934;
                c2[6, 1] := 0.1407810153019944;
                c2[6, 2] := 0.3755841316563539;
                c2[7, 1] := 0.1169627966707339;
                c2[7, 2] := 0.1937088226304455;
              elseif order == 15 then
                alpha := 0.2248854870552422;
                c1[1] := 0.4328492272335646;
                c2[1, 1] := 0.1857292591004588;
                c2[1, 2] := 0.8496337061962563;
                c2[2, 1] := 0.1808644178280136;
                c2[2, 2] := 0.8020517898136011;
                c2[3, 1] := 0.1728264404199081;
                c2[3, 2] := 0.7247449729331105;
                c2[4, 1] := 0.1616970125901954;
                c2[4, 2] := 0.6205369315943097;
                c2[5, 1] := 0.1475257264578426;
                c2[5, 2] := 0.4929612162355906;
                c2[6, 1] := 0.1301861023357119;
                c2[6, 2] := 0.3454770708040735;
                c2[7, 1] := 0.1087810777120188;
                c2[7, 2] := 0.1784526655428406;
              elseif order == 16 then
                alpha := 0.2174105053474761;
                c2[1, 1] := 0.1765637967473151;
                c2[1, 2] := 0.8377453068635511;
                c2[2, 1] := 0.1738525357503125;
                c2[2, 2] := 0.8102988957433199;
                c2[3, 1] := 0.1684627004613343;
                c2[3, 2] := 0.7563265923413258;
                c2[4, 1] := 0.1604519074815815;
                c2[4, 2] := 0.6776082294687619;
                c2[5, 1] := 0.1498828607802206;
                c2[5, 2] := 0.5766417034027680;
                c2[6, 1] := 0.1367764717792823;
                c2[6, 2] := 0.4563528264410489;
                c2[7, 1] := 0.1209810465419295;
                c2[7, 2] := 0.3193782657322374;
                c2[8, 1] := 0.1016312648007554;
                c2[8, 2] := 0.1652419227369036;
              elseif order == 17 then
                alpha := 0.2106355148193306;
                c1[1] := 0.4093223608497299;
                c2[1, 1] := 0.1664014345826274;
                c2[1, 2] := 0.8067173752345952;
                c2[2, 1] := 0.1629839591538256;
                c2[2, 2] := 0.7712924931447541;
                c2[3, 1] := 0.1573277802512491;
                c2[3, 2] := 0.7134213666303411;
                c2[4, 1] := 0.1494828185148637;
                c2[4, 2] := 0.6347841731714884;
                c2[5, 1] := 0.1394948812681826;
                c2[5, 2] := 0.5375594414619047;
                c2[6, 1] := 0.1273627583380806;
                c2[6, 2] := 0.4241608926375478;
                c2[7, 1] := 0.1129187258461290;
                c2[7, 2] := 0.2965752009703245;
                c2[8, 1] := 0.9533357359908857e-1;
                c2[8, 2] := 0.1537041700889585;
              elseif order == 18 then
                alpha := 0.2044575288651841;
                c2[1, 1] := 0.1588768571976356;
                c2[1, 2] := 0.7951914263212913;
                c2[2, 1] := 0.1569357024981854;
                c2[2, 2] := 0.7744529690772538;
                c2[3, 1] := 0.1530722206358810;
                c2[3, 2] := 0.7335304425992080;
                c2[4, 1] := 0.1473206710524167;
                c2[4, 2] := 0.6735038935387268;
                c2[5, 1] := 0.1397225420331520;
                c2[5, 2] := 0.5959151542621590;
                c2[6, 1] := 0.1303092459809849;
                c2[6, 2] := 0.5026483447894845;
                c2[7, 1] := 0.1190627367060072;
                c2[7, 2] := 0.3956893824587150;
                c2[8, 1] := 0.1058058030798994;
                c2[8, 2] := 0.2765091830730650;
                c2[9, 1] := 0.8974708108800873e-1;
                c2[9, 2] := 0.1435505288284833;
              elseif order == 19 then
                alpha := 0.1987936248083529;
                c1[1] := 0.3892259966869526;
                c2[1, 1] := 0.1506640012172225;
                c2[1, 2] := 0.7693121733774260;
                c2[2, 1] := 0.1481728062796673;
                c2[2, 2] := 0.7421133586741549;
                c2[3, 1] := 0.1440444668388838;
                c2[3, 2] := 0.6975075386214800;
                c2[4, 1] := 0.1383101628540374;
                c2[4, 2] := 0.6365464378910025;
                c2[5, 1] := 0.1310032283190998;
                c2[5, 2] := 0.5606211948462122;
                c2[6, 1] := 0.1221431166405330;
                c2[6, 2] := 0.4713530424221445;
                c2[7, 1] := 0.1116991161103884;
                c2[7, 2] := 0.3703717538617073;
                c2[8, 1] := 0.9948917351196349e-1;
                c2[8, 2] := 0.2587371155559744;
                c2[9, 1] := 0.8475989238107367e-1;
                c2[9, 2] := 0.1345537894555993;
              elseif order == 20 then
                alpha := 0.1935761760416219;
                c2[1, 1] := 0.1443871348337404;
                c2[1, 2] := 0.7584165598446141;
                c2[2, 1] := 0.1429501891353184;
                c2[2, 2] := 0.7423000962318863;
                c2[3, 1] := 0.1400877384920004;
                c2[3, 2] := 0.7104185332215555;
                c2[4, 1] := 0.1358210369491446;
                c2[4, 2] := 0.6634599783272630;
                c2[5, 1] := 0.1301773703034290;
                c2[5, 2] := 0.6024175491895959;
                c2[6, 1] := 0.1231826501439148;
                c2[6, 2] := 0.5285332736326852;
                c2[7, 1] := 0.1148465498575254;
                c2[7, 2] := 0.4431977385498628;
                c2[8, 1] := 0.1051289462376788;
                c2[8, 2] := 0.3477444062821162;
                c2[9, 1] := 0.9384622797485121e-1;
                c2[9, 2] := 0.2429038300327729;
                c2[10, 1] := 0.8028211612831444e-1;
                c2[10, 2] := 0.1265329974009533;
              elseif order == 21 then
                alpha := 0.1887494014766075;
                c1[1] := 0.3718070668941645;
                c2[1, 1] := 0.1376151928386445;
                c2[1, 2] := 0.7364290859445481;
                c2[2, 1] := 0.1357438914390695;
                c2[2, 2] := 0.7150167318935022;
                c2[3, 1] := 0.1326398453462415;
                c2[3, 2] := 0.6798001808470175;
                c2[4, 1] := 0.1283231214897678;
                c2[4, 2] := 0.6314663440439816;
                c2[5, 1] := 0.1228169159777534;
                c2[5, 2] := 0.5709353626166905;
                c2[6, 1] := 0.1161406100773184;
                c2[6, 2] := 0.4993087153571335;
                c2[7, 1] := 0.1082959649233524;
                c2[7, 2] := 0.4177766148584385;
                c2[8, 1] := 0.9923596957485723e-1;
                c2[8, 2] := 0.3274257287232124;
                c2[9, 1] := 0.8877776108724853e-1;
                c2[9, 2] := 0.2287218166767916;
                c2[10, 1] := 0.7624076527736326e-1;
                c2[10, 2] := 0.1193423971506988;
              elseif order == 22 then
                alpha := 0.1842668221199706;
                c2[1, 1] := 0.1323053462701543;
                c2[1, 2] := 0.7262446126765204;
                c2[2, 1] := 0.1312121721769772;
                c2[2, 2] := 0.7134286088450949;
                c2[3, 1] := 0.1290330911166814;
                c2[3, 2] := 0.6880287870435514;
                c2[4, 1] := 0.1257817990372067;
                c2[4, 2] := 0.6505015800059301;
                c2[5, 1] := 0.1214765261983008;
                c2[5, 2] := 0.6015107185211451;
                c2[6, 1] := 0.1161365140967959;
                c2[6, 2] := 0.5418983553698413;
                c2[7, 1] := 0.1097755171533100;
                c2[7, 2] := 0.4726370779831614;
                c2[8, 1] := 0.1023889478519956;
                c2[8, 2] := 0.3947439506537486;
                c2[9, 1] := 0.9392485861253800e-1;
                c2[9, 2] := 0.3090996703083202;
                c2[10, 1] := 0.8420273775456455e-1;
                c2[10, 2] := 0.2159561978556017;
                c2[11, 1] := 0.7257600023938262e-1;
                c2[11, 2] := 0.1128633732721116;
              elseif order == 23 then
                alpha := 0.1800893554453722;
                c1[1] := 0.3565232673929280;
                c2[1, 1] := 0.1266275171652706;
                c2[1, 2] := 0.7072778066734162;
                c2[2, 1] := 0.1251865227648538;
                c2[2, 2] := 0.6900676345785905;
                c2[3, 1] := 0.1227944815236645;
                c2[3, 2] := 0.6617011100576023;
                c2[4, 1] := 0.1194647013077667;
                c2[4, 2] := 0.6226432315773119;
                c2[5, 1] := 0.1152132989252356;
                c2[5, 2] := 0.5735222810625359;
                c2[6, 1] := 0.1100558598478487;
                c2[6, 2] := 0.5151027978024605;
                c2[7, 1] := 0.1040013558214886;
                c2[7, 2] := 0.4482410942032739;
                c2[8, 1] := 0.9704014176512626e-1;
                c2[8, 2] := 0.3738049984631116;
                c2[9, 1] := 0.8911683905758054e-1;
                c2[9, 2] := 0.2925028692588410;
                c2[10, 1] := 0.8005438265072295e-1;
                c2[10, 2] := 0.2044134600278901;
                c2[11, 1] := 0.6923832296800832e-1;
                c2[11, 2] := 0.1069984887283394;
              elseif order == 24 then
                alpha := 0.1761838665838427;
                c2[1, 1] := 0.1220804912720132;
                c2[1, 2] := 0.6978026874156063;
                c2[2, 1] := 0.1212296762358897;
                c2[2, 2] := 0.6874139794926736;
                c2[3, 1] := 0.1195328372961027;
                c2[3, 2] := 0.6667954259551859;
                c2[4, 1] := 0.1169990987333593;
                c2[4, 2] := 0.6362602049901176;
                c2[5, 1] := 0.1136409040480130;
                c2[5, 2] := 0.5962662188435553;
                c2[6, 1] := 0.1094722001757955;
                c2[6, 2] := 0.5474001634109253;
                c2[7, 1] := 0.1045052832229087;
                c2[7, 2] := 0.4903523180249535;
                c2[8, 1] := 0.9874509806025907e-1;
                c2[8, 2] := 0.4258751523524645;
                c2[9, 1] := 0.9217799943472177e-1;
                c2[9, 2] := 0.3547079765396403;
                c2[10, 1] := 0.8474633796250476e-1;
                c2[10, 2] := 0.2774145482392767;
                c2[11, 1] := 0.7627722381240495e-1;
                c2[11, 2] := 0.1939329108084139;
                c2[12, 1] := 0.6618645465422745e-1;
                c2[12, 2] := 0.1016670147947242;
              elseif order == 25 then
                alpha := 0.1725220521949266;
                c1[1] := 0.3429735385896000;
                c2[1, 1] := 0.1172525033170618;
                c2[1, 2] := 0.6812327932576614;
                c2[2, 1] := 0.1161194585333535;
                c2[2, 2] := 0.6671566071153211;
                c2[3, 1] := 0.1142375145794466;
                c2[3, 2] := 0.6439167855053158;
                c2[4, 1] := 0.1116157454252308;
                c2[4, 2] := 0.6118378416180135;
                c2[5, 1] := 0.1082654809459177;
                c2[5, 2] := 0.5713609763370088;
                c2[6, 1] := 0.1041985674230918;
                c2[6, 2] := 0.5230289949762722;
                c2[7, 1] := 0.9942439308123559e-1;
                c2[7, 2] := 0.4674627926041906;
                c2[8, 1] := 0.9394453593830893e-1;
                c2[8, 2] := 0.4053226688298811;
                c2[9, 1] := 0.8774221237222533e-1;
                c2[9, 2] := 0.3372372276379071;
                c2[10, 1] := 0.8075839512216483e-1;
                c2[10, 2] := 0.2636485508005428;
                c2[11, 1] := 0.7282483286646764e-1;
                c2[11, 2] := 0.1843801345273085;
                c2[12, 1] := 0.6338571166846652e-1;
                c2[12, 2] := 0.9680153764737715e-1;
              elseif order == 26 then
                alpha := 0.1690795702796737;
                c2[1, 1] := 0.1133168695796030;
                c2[1, 2] := 0.6724297955493932;
                c2[2, 1] := 0.1126417845769961;
                c2[2, 2] := 0.6638709519790540;
                c2[3, 1] := 0.1112948749545606;
                c2[3, 2] := 0.6468652038763624;
                c2[4, 1] := 0.1092823986944244;
                c2[4, 2] := 0.6216337070799265;
                c2[5, 1] := 0.1066130386697976;
                c2[5, 2] := 0.5885011413992190;
                c2[6, 1] := 0.1032969057045413;
                c2[6, 2] := 0.5478864278297548;
                c2[7, 1] := 0.9934388184210715e-1;
                c2[7, 2] := 0.5002885306054287;
                c2[8, 1] := 0.9476081523436283e-1;
                c2[8, 2] := 0.4462644847551711;
                c2[9, 1] := 0.8954648464575577e-1;
                c2[9, 2] := 0.3863930785049522;
                c2[10, 1] := 0.8368166847159917e-1;
                c2[10, 2] := 0.3212074592527143;
                c2[11, 1] := 0.7710664731701103e-1;
                c2[11, 2] := 0.2510470347119383;
                c2[12, 1] := 0.6965807988411425e-1;
                c2[12, 2] := 0.1756419294111342;
                c2[13, 1] := 0.6080674930548766e-1;
                c2[13, 2] := 0.9234535279274277e-1;
              elseif order == 27 then
                alpha := 0.1658353543067995;
                c1[1] := 0.3308543720638957;
                c2[1, 1] := 0.1091618578712746;
                c2[1, 2] := 0.6577977071169651;
                c2[2, 1] := 0.1082549561495043;
                c2[2, 2] := 0.6461121666520275;
                c2[3, 1] := 0.1067479247890451;
                c2[3, 2] := 0.6267937760991321;
                c2[4, 1] := 0.1046471079537577;
                c2[4, 2] := 0.6000750116745808;
                c2[5, 1] := 0.1019605976654259;
                c2[5, 2] := 0.5662734183049320;
                c2[6, 1] := 0.9869726954433709e-1;
                c2[6, 2] := 0.5257827234948534;
                c2[7, 1] := 0.9486520934132483e-1;
                c2[7, 2] := 0.4790595019077763;
                c2[8, 1] := 0.9046906518775348e-1;
                c2[8, 2] := 0.4266025862147336;
                c2[9, 1] := 0.8550529998276152e-1;
                c2[9, 2] := 0.3689188223512328;
                c2[10, 1] := 0.7995282239306020e-1;
                c2[10, 2] := 0.3064589322702932;
                c2[11, 1] := 0.7375174596252882e-1;
                c2[11, 2] := 0.2394754504667310;
                c2[12, 1] := 0.6674377263329041e-1;
                c2[12, 2] := 0.1676223546666024;
                c2[13, 1] := 0.5842458027529246e-1;
                c2[13, 2] := 0.8825044329219431e-1;
              elseif order == 28 then
                alpha := 0.1627710671942929;
                c2[1, 1] := 0.1057232656113488;
                c2[1, 2] := 0.6496161226860832;
                c2[2, 1] := 0.1051786825724864;
                c2[2, 2] := 0.6424661279909941;
                c2[3, 1] := 0.1040917964935006;
                c2[3, 2] := 0.6282470268918791;
                c2[4, 1] := 0.1024670101953951;
                c2[4, 2] := 0.6071189030701136;
                c2[5, 1] := 0.1003105109519892;
                c2[5, 2] := 0.5793175191747016;
                c2[6, 1] := 0.9762969425430802e-1;
                c2[6, 2] := 0.5451486608855443;
                c2[7, 1] := 0.9443223803058400e-1;
                c2[7, 2] := 0.5049796971628137;
                c2[8, 1] := 0.9072460982036488e-1;
                c2[8, 2] := 0.4592270546572523;
                c2[9, 1] := 0.8650956423253280e-1;
                c2[9, 2] := 0.4083368605952977;
                c2[10, 1] := 0.8178165740374893e-1;
                c2[10, 2] := 0.3527525188880655;
                c2[11, 1] := 0.7651838885868020e-1;
                c2[11, 2] := 0.2928534570013572;
                c2[12, 1] := 0.7066010532447490e-1;
                c2[12, 2] := 0.2288185204390681;
                c2[13, 1] := 0.6405358596145789e-1;
                c2[13, 2] := 0.1602396172588190;
                c2[14, 1] := 0.5621780070227172e-1;
                c2[14, 2] := 0.8447589564915071e-1;
              elseif order == 29 then
                alpha := 0.1598706626277596;
                c1[1] := 0.3199314513011623;
                c2[1, 1] := 0.1021101032532951;
                c2[1, 2] := 0.6365758882240111;
                c2[2, 1] := 0.1013729819392774;
                c2[2, 2] := 0.6267495975736321;
                c2[3, 1] := 0.1001476175660628;
                c2[3, 2] := 0.6104876178266819;
                c2[4, 1] := 0.9843854640428316e-1;
                c2[4, 2] := 0.5879603139195113;
                c2[5, 1] := 0.9625164534591696e-1;
                c2[5, 2] := 0.5594012291050210;
                c2[6, 1] := 0.9359356960417668e-1;
                c2[6, 2] := 0.5251016150410664;
                c2[7, 1] := 0.9047086748649986e-1;
                c2[7, 2] := 0.4854024475590397;
                c2[8, 1] := 0.8688856407189167e-1;
                c2[8, 2] := 0.4406826457109709;
                c2[9, 1] := 0.8284779224069856e-1;
                c2[9, 2] := 0.3913408089298914;
                c2[10, 1] := 0.7834154620997181e-1;
                c2[10, 2] := 0.3377643999400627;
                c2[11, 1] := 0.7334628941928766e-1;
                c2[11, 2] := 0.2802710651919946;
                c2[12, 1] := 0.6780290487362146e-1;
                c2[12, 2] := 0.2189770008083379;
                c2[13, 1] := 0.6156321231528423e-1;
                c2[13, 2] := 0.1534235999306070;
                c2[14, 1] := 0.5416797446761512e-1;
                c2[14, 2] := 0.8098664736760292e-1;
              elseif order == 30 then
                alpha := 0.1571200296252450;
                c2[1, 1] := 0.9908074847842124e-1;
                c2[1, 2] := 0.6289618807831557;
                c2[2, 1] := 0.9863509708328196e-1;
                c2[2, 2] := 0.6229164525571278;
                c2[3, 1] := 0.9774542692037148e-1;
                c2[3, 2] := 0.6108853364240036;
                c2[4, 1] := 0.9641490581986484e-1;
                c2[4, 2] := 0.5929869253412513;
                c2[5, 1] := 0.9464802912225441e-1;
                c2[5, 2] := 0.5693960175547550;
                c2[6, 1] := 0.9245027206218041e-1;
                c2[6, 2] := 0.5403402396359503;
                c2[7, 1] := 0.8982754584112941e-1;
                c2[7, 2] := 0.5060948065875106;
                c2[8, 1] := 0.8678535291732599e-1;
                c2[8, 2] := 0.4669749797983789;
                c2[9, 1] := 0.8332744242052199e-1;
                c2[9, 2] := 0.4233249626334694;
                c2[10, 1] := 0.7945356393775309e-1;
                c2[10, 2] := 0.3755006094498054;
                c2[11, 1] := 0.7515543969833788e-1;
                c2[11, 2] := 0.3238400339292700;
                c2[12, 1] := 0.7040879901685638e-1;
                c2[12, 2] := 0.2686072427439079;
                c2[13, 1] := 0.6515528854010540e-1;
                c2[13, 2] := 0.2098650589782619;
                c2[14, 1] := 0.5925168237177876e-1;
                c2[14, 2] := 0.1471138832654873;
                c2[15, 1] := 0.5225913954211672e-1;
                c2[15, 2] := 0.7775248839507864e-1;
              elseif order == 31 then
                alpha := 0.1545067022920929;
                c1[1] := 0.3100206996451866;
                c2[1, 1] := 0.9591020358831668e-1;
                c2[1, 2] := 0.6172474793293396;
                c2[2, 1] := 0.9530301275601203e-1;
                c2[2, 2] := 0.6088916323460413;
                c2[3, 1] := 0.9429332655402368e-1;
                c2[3, 2] := 0.5950511595503025;
                c2[4, 1] := 0.9288445429894548e-1;
                c2[4, 2] := 0.5758534119053522;
                c2[5, 1] := 0.9108073420087422e-1;
                c2[5, 2] := 0.5514734636081183;
                c2[6, 1] := 0.8888719137536870e-1;
                c2[6, 2] := 0.5221306199481831;
                c2[7, 1] := 0.8630901440239650e-1;
                c2[7, 2] := 0.4880834248148061;
                c2[8, 1] := 0.8335074993373294e-1;
                c2[8, 2] := 0.4496225358496770;
                c2[9, 1] := 0.8001502494376102e-1;
                c2[9, 2] := 0.4070602306679052;
                c2[10, 1] := 0.7630041338037624e-1;
                c2[10, 2] := 0.3607139804818122;
                c2[11, 1] := 0.7219760885744920e-1;
                c2[11, 2] := 0.3108783301229550;
                c2[12, 1] := 0.6768185077153345e-1;
                c2[12, 2] := 0.2577706252514497;
                c2[13, 1] := 0.6269571766328638e-1;
                c2[13, 2] := 0.2014081375889921;
                c2[14, 1] := 0.5710081766945065e-1;
                c2[14, 2] := 0.1412581515841926;
                c2[15, 1] := 0.5047740914807019e-1;
                c2[15, 2] := 0.7474725873250158e-1;
              elseif order == 32 then
                alpha := 0.1520196210848210;
                c2[1, 1] := 0.9322163554339406e-1;
                c2[1, 2] := 0.6101488690506050;
                c2[2, 1] := 0.9285233997694042e-1;
                c2[2, 2] := 0.6049832320721264;
                c2[3, 1] := 0.9211494244473163e-1;
                c2[3, 2] := 0.5946969295569034;
                c2[4, 1] := 0.9101176786042449e-1;
                c2[4, 2] := 0.5793791854364477;
                c2[5, 1] := 0.8954614071360517e-1;
                c2[5, 2] := 0.5591619969234026;
                c2[6, 1] := 0.8772216763680164e-1;
                c2[6, 2] := 0.5342177994699602;
                c2[7, 1] := 0.8554440426912734e-1;
                c2[7, 2] := 0.5047560942986598;
                c2[8, 1] := 0.8301735302045588e-1;
                c2[8, 2] := 0.4710187048140929;
                c2[9, 1] := 0.8014469519188161e-1;
                c2[9, 2] := 0.4332730387207936;
                c2[10, 1] := 0.7692807528893225e-1;
                c2[10, 2] := 0.3918021436411035;
                c2[11, 1] := 0.7336507157284898e-1;
                c2[11, 2] := 0.3468890521471250;
                c2[12, 1] := 0.6944555312763458e-1;
                c2[12, 2] := 0.2987898029050460;
                c2[13, 1] := 0.6514446669420571e-1;
                c2[13, 2] := 0.2476810747407199;
                c2[14, 1] := 0.6040544477732702e-1;
                c2[14, 2] := 0.1935412053397663;
                c2[15, 1] := 0.5509478650672775e-1;
                c2[15, 2] := 0.1358108994174911;
                c2[16, 1] := 0.4881064725720192e-1;
                c2[16, 2] := 0.7194819894416505e-1;
              elseif order == 33 then
                alpha := 0.1496489351138032;
                c1[1] := 0.3009752799176432;
                c2[1, 1] := 0.9041725460994505e-1;
                c2[1, 2] := 0.5995521047364046;
                c2[2, 1] := 0.8991117804113002e-1;
                c2[2, 2] := 0.5923764112099496;
                c2[3, 1] := 0.8906941547422532e-1;
                c2[3, 2] := 0.5804822013853129;
                c2[4, 1] := 0.8789442491445575e-1;
                c2[4, 2] := 0.5639663528946501;
                c2[5, 1] := 0.8638945831033775e-1;
                c2[5, 2] := 0.5429623519607796;
                c2[6, 1] := 0.8455834602616358e-1;
                c2[6, 2] := 0.5176379938389326;
                c2[7, 1] := 0.8240517431382334e-1;
                c2[7, 2] := 0.4881921474066189;
                c2[8, 1] := 0.7993380417355076e-1;
                c2[8, 2] := 0.4548502528082586;
                c2[9, 1] := 0.7714713890732801e-1;
                c2[9, 2] := 0.4178579388038483;
                c2[10, 1] := 0.7404596598181127e-1;
                c2[10, 2] := 0.3774715722484659;
                c2[11, 1] := 0.7062702339160462e-1;
                c2[11, 2] := 0.3339432938810453;
                c2[12, 1] := 0.6687952672391507e-1;
                c2[12, 2] := 0.2874950693388235;
                c2[13, 1] := 0.6277828912909767e-1;
                c2[13, 2] := 0.2382680702894708;
                c2[14, 1] := 0.5826808305383988e-1;
                c2[14, 2] := 0.1862073169968455;
                c2[15, 1] := 0.5321974125363517e-1;
                c2[15, 2] := 0.1307323751236313;
                c2[16, 1] := 0.4724820282032780e-1;
                c2[16, 2] := 0.6933542082177094e-1;
              elseif order == 34 then
                alpha := 0.1473858373968463;
                c2[1, 1] := 0.8801537152275983e-1;
                c2[1, 2] := 0.5929204288972172;
                c2[2, 1] := 0.8770594341007476e-1;
                c2[2, 2] := 0.5884653382247518;
                c2[3, 1] := 0.8708797598072095e-1;
                c2[3, 2] := 0.5795895850253119;
                c2[4, 1] := 0.8616320590689187e-1;
                c2[4, 2] := 0.5663615383647170;
                c2[5, 1] := 0.8493413175570858e-1;
                c2[5, 2] := 0.5488825092350877;
                c2[6, 1] := 0.8340387368687513e-1;
                c2[6, 2] := 0.5272851839324592;
                c2[7, 1] := 0.8157596213131521e-1;
                c2[7, 2] := 0.5017313864372913;
                c2[8, 1] := 0.7945402670834270e-1;
                c2[8, 2] := 0.4724089864574216;
                c2[9, 1] := 0.7704133559556429e-1;
                c2[9, 2] := 0.4395276256463053;
                c2[10, 1] := 0.7434009635219704e-1;
                c2[10, 2] := 0.4033126590648964;
                c2[11, 1] := 0.7135035113853376e-1;
                c2[11, 2] := 0.3639961488919042;
                c2[12, 1] := 0.6806813160738834e-1;
                c2[12, 2] := 0.3218025212900124;
                c2[13, 1] := 0.6448214312000864e-1;
                c2[13, 2] := 0.2769235521088158;
                c2[14, 1] := 0.6056719318430530e-1;
                c2[14, 2] := 0.2294693573271038;
                c2[15, 1] := 0.5626925196925040e-1;
                c2[15, 2] := 0.1793564218840015;
                c2[16, 1] := 0.5146352031547277e-1;
                c2[16, 2] := 0.1259877129326412;
                c2[17, 1] := 0.4578069074410591e-1;
                c2[17, 2] := 0.6689147319568768e-1;
              elseif order == 35 then
                alpha := 0.1452224267615486;
                c1[1] := 0.2926764667564367;
                c2[1, 1] := 0.8551731299267280e-1;
                c2[1, 2] := 0.5832758214629523;
                c2[2, 1] := 0.8509109732853060e-1;
                c2[2, 2] := 0.5770596582643844;
                c2[3, 1] := 0.8438201446671953e-1;
                c2[3, 2] := 0.5667497616665494;
                c2[4, 1] := 0.8339191981579831e-1;
                c2[4, 2] := 0.5524209816238369;
                c2[5, 1] := 0.8212328610083385e-1;
                c2[5, 2] := 0.5341766459916322;
                c2[6, 1] := 0.8057906332198853e-1;
                c2[6, 2] := 0.5121470053512750;
                c2[7, 1] := 0.7876247299954955e-1;
                c2[7, 2] := 0.4864870722254752;
                c2[8, 1] := 0.7667670879950268e-1;
                c2[8, 2] := 0.4573736721705665;
                c2[9, 1] := 0.7432449556218945e-1;
                c2[9, 2] := 0.4250013835198991;
                c2[10, 1] := 0.7170742126011575e-1;
                c2[10, 2] := 0.3895767735915445;
                c2[11, 1] := 0.6882488171701314e-1;
                c2[11, 2] := 0.3513097926737368;
                c2[12, 1] := 0.6567231746957568e-1;
                c2[12, 2] := 0.3103999917596611;
                c2[13, 1] := 0.6223804362223595e-1;
                c2[13, 2] := 0.2670123611280899;
                c2[14, 1] := 0.5849696460782910e-1;
                c2[14, 2] := 0.2212298104867592;
                c2[15, 1] := 0.5439628409499822e-1;
                c2[15, 2] := 0.1729443731341637;
                c2[16, 1] := 0.4981540179136920e-1;
                c2[16, 2] := 0.1215462157134930;
                c2[17, 1] := 0.4439981033536435e-1;
                c2[17, 2] := 0.6460098363520967e-1;
              elseif order == 36 then
                alpha := 0.1431515914458580;
                c2[1, 1] := 0.8335881847130301e-1;
                c2[1, 2] := 0.5770670512160201;
                c2[2, 1] := 0.8309698922852212e-1;
                c2[2, 2] := 0.5731929100172432;
                c2[3, 1] := 0.8257400347039723e-1;
                c2[3, 2] := 0.5654713811993058;
                c2[4, 1] := 0.8179117911600136e-1;
                c2[4, 2] := 0.5539556343603020;
                c2[5, 1] := 0.8075042173126963e-1;
                c2[5, 2] := 0.5387245649546684;
                c2[6, 1] := 0.7945413151258206e-1;
                c2[6, 2] := 0.5198817177723069;
                c2[7, 1] := 0.7790506514288866e-1;
                c2[7, 2] := 0.4975537629595409;
                c2[8, 1] := 0.7610613635339480e-1;
                c2[8, 2] := 0.4718884193866789;
                c2[9, 1] := 0.7406012816626425e-1;
                c2[9, 2] := 0.4430516443136726;
                c2[10, 1] := 0.7176927060205631e-1;
                c2[10, 2] := 0.4112237708115829;
                c2[11, 1] := 0.6923460172504251e-1;
                c2[11, 2] := 0.3765940116389730;
                c2[12, 1] := 0.6645495833489556e-1;
                c2[12, 2] := 0.3393522147815403;
                c2[13, 1] := 0.6342528888937094e-1;
                c2[13, 2] := 0.2996755899575573;
                c2[14, 1] := 0.6013361864949449e-1;
                c2[14, 2] := 0.2577053294053830;
                c2[15, 1] := 0.5655503081322404e-1;
                c2[15, 2] := 0.2135004731531631;
                c2[16, 1] := 0.5263798119559069e-1;
                c2[16, 2] := 0.1669320999865636;
                c2[17, 1] := 0.4826589873626196e-1;
                c2[17, 2] := 0.1173807590715484;
                c2[18, 1] := 0.4309819397289806e-1;
                c2[18, 2] := 0.6245036108880222e-1;
              elseif order == 37 then
                alpha := 0.1411669104782917;
                c1[1] := 0.2850271036215707;
                c2[1, 1] := 0.8111958235023328e-1;
                c2[1, 2] := 0.5682412610563970;
                c2[2, 1] := 0.8075727567979578e-1;
                c2[2, 2] := 0.5628142923227016;
                c2[3, 1] := 0.8015440554413301e-1;
                c2[3, 2] := 0.5538087696879930;
                c2[4, 1] := 0.7931239302677386e-1;
                c2[4, 2] := 0.5412833323304460;
                c2[5, 1] := 0.7823314328639347e-1;
                c2[5, 2] := 0.5253190555393968;
                c2[6, 1] := 0.7691895211595101e-1;
                c2[6, 2] := 0.5060183741977191;
                c2[7, 1] := 0.7537237072011853e-1;
                c2[7, 2] := 0.4835036020049034;
                c2[8, 1] := 0.7359601294804538e-1;
                c2[8, 2] := 0.4579149413954837;
                c2[9, 1] := 0.7159227884849299e-1;
                c2[9, 2] := 0.4294078049978829;
                c2[10, 1] := 0.6936295002846032e-1;
                c2[10, 2] := 0.3981491350382047;
                c2[11, 1] := 0.6690857785828917e-1;
                c2[11, 2] := 0.3643121502867948;
                c2[12, 1] := 0.6422751692085542e-1;
                c2[12, 2] := 0.3280684291406284;
                c2[13, 1] := 0.6131430866206096e-1;
                c2[13, 2] := 0.2895750997170303;
                c2[14, 1] := 0.5815677249570920e-1;
                c2[14, 2] := 0.2489521814805720;
                c2[15, 1] := 0.5473023527947980e-1;
                c2[15, 2] := 0.2062377435955363;
                c2[16, 1] := 0.5098441033167034e-1;
                c2[16, 2] := 0.1612849131645336;
                c2[17, 1] := 0.4680658811093562e-1;
                c2[17, 2] := 0.1134672937045305;
                c2[18, 1] := 0.4186928031694695e-1;
                c2[18, 2] := 0.6042754777339966e-1;
              elseif order == 38 then
                alpha := 0.1392625697140030;
                c2[1, 1] := 0.7916943373658329e-1;
                c2[1, 2] := 0.5624158631591745;
                c2[2, 1] := 0.7894592250257840e-1;
                c2[2, 2] := 0.5590219398777304;
                c2[3, 1] := 0.7849941672384930e-1;
                c2[3, 2] := 0.5522551628416841;
                c2[4, 1] := 0.7783093084875645e-1;
                c2[4, 2] := 0.5421574325808380;
                c2[5, 1] := 0.7694193770482690e-1;
                c2[5, 2] := 0.5287909941093643;
                c2[6, 1] := 0.7583430534712885e-1;
                c2[6, 2] := 0.5122376814029880;
                c2[7, 1] := 0.7451020436122948e-1;
                c2[7, 2] := 0.4925978555548549;
                c2[8, 1] := 0.7297197617673508e-1;
                c2[8, 2] := 0.4699889739625235;
                c2[9, 1] := 0.7122194706992953e-1;
                c2[9, 2] := 0.4445436860615774;
                c2[10, 1] := 0.6926216260386816e-1;
                c2[10, 2] := 0.4164072786327193;
                c2[11, 1] := 0.6709399961255503e-1;
                c2[11, 2] := 0.3857341621868851;
                c2[12, 1] := 0.6471757977022456e-1;
                c2[12, 2] := 0.3526828388476838;
                c2[13, 1] := 0.6213084287116965e-1;
                c2[13, 2] := 0.3174082831364342;
                c2[14, 1] := 0.5932799638550641e-1;
                c2[14, 2] := 0.2800495563550299;
                c2[15, 1] := 0.5629672408524944e-1;
                c2[15, 2] := 0.2407078154782509;
                c2[16, 1] := 0.5301264751544952e-1;
                c2[16, 2] := 0.1994026830553859;
                c2[17, 1] := 0.4942673259817896e-1;
                c2[17, 2] := 0.1559719194038917;
                c2[18, 1] := 0.4542996716979947e-1;
                c2[18, 2] := 0.1097844277878470;
                c2[19, 1] := 0.4070720755433961e-1;
                c2[19, 2] := 0.5852181110523043e-1;
              elseif order == 39 then
                alpha := 0.1374332900196804;
                c1[1] := 0.2779468246419593;
                c2[1, 1] := 0.7715084161825772e-1;
                c2[1, 2] := 0.5543001331300056;
                c2[2, 1] := 0.7684028301163326e-1;
                c2[2, 2] := 0.5495289890712267;
                c2[3, 1] := 0.7632343924866024e-1;
                c2[3, 2] := 0.5416083298429741;
                c2[4, 1] := 0.7560141319808483e-1;
                c2[4, 2] := 0.5305846713929198;
                c2[5, 1] := 0.7467569064745969e-1;
                c2[5, 2] := 0.5165224112570647;
                c2[6, 1] := 0.7354807648551346e-1;
                c2[6, 2] := 0.4995030679271456;
                c2[7, 1] := 0.7222060351121389e-1;
                c2[7, 2] := 0.4796242430956156;
                c2[8, 1] := 0.7069540462458585e-1;
                c2[8, 2] := 0.4569982440368368;
                c2[9, 1] := 0.6897453353492381e-1;
                c2[9, 2] := 0.4317502624832354;
                c2[10, 1] := 0.6705970959388781e-1;
                c2[10, 2] := 0.4040159353969854;
                c2[11, 1] := 0.6495194541066725e-1;
                c2[11, 2] := 0.3739379843169939;
                c2[12, 1] := 0.6265098412417610e-1;
                c2[12, 2] := 0.3416613843816217;
                c2[13, 1] := 0.6015440984955930e-1;
                c2[13, 2] := 0.3073260166338746;
                c2[14, 1] := 0.5745615876877304e-1;
                c2[14, 2] := 0.2710546723961181;
                c2[15, 1] := 0.5454383762391338e-1;
                c2[15, 2] := 0.2329316824061170;
                c2[16, 1] := 0.5139340231935751e-1;
                c2[16, 2] := 0.1929604256043231;
                c2[17, 1] := 0.4795705862458131e-1;
                c2[17, 2] := 0.1509655259246037;
                c2[18, 1] := 0.4412933231935506e-1;
                c2[18, 2] := 0.1063130748962878;
                c2[19, 1] := 0.3960672309405603e-1;
                c2[19, 2] := 0.5672356837211527e-1;
              elseif order == 40 then
                alpha := 0.1356742655825434;
                c2[1, 1] := 0.7538038374294594e-1;
                c2[1, 2] := 0.5488228264329617;
                c2[2, 1] := 0.7518806529402738e-1;
                c2[2, 2] := 0.5458297722483311;
                c2[3, 1] := 0.7480383050347119e-1;
                c2[3, 2] := 0.5398604576730540;
                c2[4, 1] := 0.7422847031965465e-1;
                c2[4, 2] := 0.5309482987446206;
                c2[5, 1] := 0.7346313704205006e-1;
                c2[5, 2] := 0.5191429845322307;
                c2[6, 1] := 0.7250930053201402e-1;
                c2[6, 2] := 0.5045099368431007;
                c2[7, 1] := 0.7136868456879621e-1;
                c2[7, 2] := 0.4871295553902607;
                c2[8, 1] := 0.7004317764946634e-1;
                c2[8, 2] := 0.4670962098860498;
                c2[9, 1] := 0.6853470921527828e-1;
                c2[9, 2] := 0.4445169164956202;
                c2[10, 1] := 0.6684507689945471e-1;
                c2[10, 2] := 0.4195095960479698;
                c2[11, 1] := 0.6497570123412630e-1;
                c2[11, 2] := 0.3922007419030645;
                c2[12, 1] := 0.6292726794917847e-1;
                c2[12, 2] := 0.3627221993494397;
                c2[13, 1] := 0.6069918741663154e-1;
                c2[13, 2] := 0.3312065181294388;
                c2[14, 1] := 0.5828873983769410e-1;
                c2[14, 2] := 0.2977798532686911;
                c2[15, 1] := 0.5568964389813015e-1;
                c2[15, 2] := 0.2625503293999835;
                c2[16, 1] := 0.5288947816690705e-1;
                c2[16, 2] := 0.2255872486520188;
                c2[17, 1] := 0.4986456327645859e-1;
                c2[17, 2] := 0.1868796731919594;
                c2[18, 1] := 0.4656832613054458e-1;
                c2[18, 2] := 0.1462410193532463;
                c2[19, 1] := 0.4289867647614935e-1;
                c2[19, 2] := 0.1030361558710747;
                c2[20, 1] := 0.3856310684054106e-1;
                c2[20, 2] := 0.5502423832293889e-1;
              elseif order == 41 then
                alpha := 0.1339811106984253;
                c1[1] := 0.2713685065531391;
                c2[1, 1] := 0.7355140275160984e-1;
                c2[1, 2] := 0.5413274778282860;
                c2[2, 1] := 0.7328319082267173e-1;
                c2[2, 2] := 0.5371064088294270;
                c2[3, 1] := 0.7283676160772547e-1;
                c2[3, 2] := 0.5300963437270770;
                c2[4, 1] := 0.7221298133014343e-1;
                c2[4, 2] := 0.5203345998371490;
                c2[5, 1] := 0.7141302173623395e-1;
                c2[5, 2] := 0.5078728971879841;
                c2[6, 1] := 0.7043831559982149e-1;
                c2[6, 2] := 0.4927768111819803;
                c2[7, 1] := 0.6929049381827268e-1;
                c2[7, 2] := 0.4751250308594139;
                c2[8, 1] := 0.6797129849758392e-1;
                c2[8, 2] := 0.4550083840638406;
                c2[9, 1] := 0.6648246325101609e-1;
                c2[9, 2] := 0.4325285673076087;
                c2[10, 1] := 0.6482554675958526e-1;
                c2[10, 2] := 0.4077964789091151;
                c2[11, 1] := 0.6300169683004558e-1;
                c2[11, 2] := 0.3809299858742483;
                c2[12, 1] := 0.6101130648543355e-1;
                c2[12, 2] := 0.3520508315700898;
                c2[13, 1] := 0.5885349417435808e-1;
                c2[13, 2] := 0.3212801560701271;
                c2[14, 1] := 0.5652528148656809e-1;
                c2[14, 2] := 0.2887316252774887;
                c2[15, 1] := 0.5402021575818373e-1;
                c2[15, 2] := 0.2545001287790888;
                c2[16, 1] := 0.5132588802608274e-1;
                c2[16, 2] := 0.2186415296842951;
                c2[17, 1] := 0.4841900639702602e-1;
                c2[17, 2] := 0.1811322622296060;
                c2[18, 1] := 0.4525419574485134e-1;
                c2[18, 2] := 0.1417762065404688;
                c2[19, 1] := 0.4173260173087802e-1;
                c2[19, 2] := 0.9993834530966510e-1;
                c2[20, 1] := 0.3757210572966463e-1;
                c2[20, 2] := 0.5341611499960143e-1;
              else
                Streams.error("Input argument order (= " + String(order) +
                  ") of Bessel filter is not in the range 1..41");
              end if;

              annotation (Documentation(info="<html><p>The transfer function H(p) of a <em>n</em> 'th order Bessel filter is given by</p>
<blockquote><pre>
        Bn(0)
H(p) = -------
        Bn(p)
</pre></blockquote>
<p>with the denominator polynomial</p>
<blockquote><pre>
         n             n  (2n - k)!       p^k
Bn(p) = sum c_k*p^k = sum ----------- * -------   (1)
        k=0           k=0 (n - k)!k!    2^(n-k)
</pre></blockquote>
<p>and the numerator</p>
<blockquote><pre>
               (2n)!     1
Bn(0) = c_0 = ------- * ---- .                    (2)
                n!      2^n
</pre></blockquote>
<p>Although the coefficients c_k are integer numbers, it is not advisable to use the
polynomials in an unfactorized form because the coefficients are fast growing with order
n (c_0 is approximately 0.3e24 and 0.8e59 for order n=20 and order n=40
respectively).</p>

<p>Therefore, the polynomial Bn(p) is factorized to first and second order polynomials with
real coefficients corresponding to zeros and poles representation that is used in this library.</p>

<p>The function returns the coefficients which resulted from factorization of the normalized transfer function</p>
<blockquote><pre>
H'(p') = H(p),  p' = p/w0
</pre></blockquote>
<p>as well as</p>
<blockquote><pre>
alpha = 1/w0
</pre></blockquote>
<p>the reciprocal of the cut of frequency w0 where the gain of the transfer function is
decreased 3dB.</p>

<p>Both, coefficients and cut off frequency were calculated symbolically and were eventually evaluated
with high precision calculation. The results were stored in this function as real
numbers.</p>

<h4>Calculation of normalized Bessel filter coefficients</h4>
<p>Equation</p>
<blockquote><pre>
abs(H(j*w0)) = abs(Bn(0)/Bn(j*w0)) = 10^(-3/20)
</pre></blockquote>
<p>which must be fulfilled for cut off frequency w = w0 leads to</p>
<blockquote><pre>
[Re(Bn(j*w0))]^2 + [Im(Bn(j*w0))]^2 - (Bn(0)^2)*10^(3/10) = 0
</pre></blockquote>
<p>which has exactly one real solution w0 for each order n. This solutions of w0 are
calculated symbolically first and evaluated by using high precise values of the
coefficients c_k calculated by following (1) and (2).</p>

<p>With w0, the coefficients of the factorized polynomial can be computed by calculating the
zeros of the denominator polynomial</p>
<blockquote><pre>
        n
Bn(p) = sum w0^k*c_k*(p/w0)^k
        k=0
</pre></blockquote>
<p>of the normalized transfer function H'(p'). There exist n/2 of conjugate complex
pairs of zeros (beta +-j*gamma) if n is even and one additional real zero (alpha) if n is
odd. Finally, the coefficients a, b1_k, b2_k of the polynomials</p>
<blockquote><pre>
a*p + 1,  n is odd
</pre></blockquote>
<p>and</p>
<blockquote><pre>
b2_k*p^2 + b1_k*p + 1,   k = 1,... div(n,2)
</pre></blockquote>
<p>results from</p>
<blockquote><pre>
a = -1/alpha
</pre></blockquote>
<p>and</p>
<blockquote><pre>
b2_k = 1/(beta_k^2 + gamma_k^2) b1_k = -2*beta_k/(beta_k^2 + gamma_k^2)
</pre></blockquote>
</html>"));
            end BesselBaseCoefficients;

            function toHighestPowerOne
              "Transform filter to form with highest power of s equal 1"
              extends Modelica.Icons.Function;

              input Real den1[:] "[s] coefficients of polynomials (den1[i]*s + 1)";
              input Real den2[:,2]
                "[s^2, s] coefficients of polynomials (den2[i,1]*s^2 + den2[i,2]*s + 1)";
              output Real cr[size(den1, 1)]
                "[s^0] coefficients of polynomials cr[i]*(s+1/cr[i])";
              output Real c0[size(den2, 1)]
                "[s^0] coefficients of polynomials (s^2 + (den2[i,2]/den2[i,1])*s + (1/den2[i,1]))";
              output Real c1[size(den2, 1)]
                "[s^1] coefficients of polynomials (s^2 + (den2[i,2]/den2[i,1])*s + (1/den2[i,1]))";
            algorithm
              for i in 1:size(den1, 1) loop
                cr[i] := 1/den1[i];
              end for;

              for i in 1:size(den2, 1) loop
                c1[i] := den2[i, 2]/den2[i, 1];
                c0[i] := 1/den2[i, 1];
              end for;
            end toHighestPowerOne;

            function normalizationFactor
              "Compute correction factor of low pass filter such that amplitude at cut-off frequency is -3db (=10^(-3/20) = 0.70794...)"
              extends Modelica.Icons.Function;

              import Modelica.Utilities.Streams;

              input Real c1[:]
                "[p] coefficients of denominator polynomials (c1[i}*p + 1)";
              input Real c2[:,2]
                "[p^2, p] coefficients of denominator polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
              output Real alpha "Correction factor (replace p by alpha*p)";
            protected
              Real alpha_min;
              Real alpha_max;

              function normalizationResidue
                "Residue of correction factor computation"
                extends Modelica.Icons.Function;
                input Real c1[:]
                  "[p] coefficients of denominator polynomials (c1[i]*p + 1)";
                input Real c2[:,2]
                  "[p^2, p] coefficients of denominator polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
                input Real alpha;
                output Real residue;
              protected
                constant Real beta= 10^(-3/20)
                  "Amplitude of -3db required, i.e., -3db = 20*log(beta)";
                Real cc1;
                Real cc2;
                Real p;
                Real alpha2=alpha*alpha;
                Real alpha4=alpha2*alpha2;
                Real A2=1.0;
              algorithm
                assert(size(c1,1) <= 1, "Internal error 2 (should not occur)");
                if size(c1, 1) == 1 then
                  cc1 := c1[1]*c1[1];
                  p := 1 + cc1*alpha2;
                  A2 := A2*p;
                end if;
                for i in 1:size(c2, 1) loop
                  cc1 := c2[i, 2]*c2[i, 2] - 2*c2[i, 1];
                  cc2 := c2[i, 1]*c2[i, 1];
                  p := 1 + cc1*alpha2 + cc2*alpha4;
                  A2 := A2*p;
                end for;
                residue := 1/sqrt(A2) - beta;
              end normalizationResidue;

              function findInterval "Find interval for the root"
                extends Modelica.Icons.Function;
                input Real c1[:]
                  "[p] coefficients of denominator polynomials (a*p + 1)";
                input Real c2[:,2]
                  "[p^2, p] coefficients of denominator polynomials (b*p^2 + a*p + 1)";
                output Real alpha_min;
                output Real alpha_max;
              protected
                Real alpha = 1.0;
                Real residue;
              algorithm
                alpha_min :=0;
                residue := normalizationResidue(c1, c2, alpha);
                if residue < 0 then
                   alpha_max :=alpha;
                else
                   while residue >= 0 loop
                      alpha := 1.1*alpha;
                      residue := normalizationResidue(c1, c2, alpha);
                   end while;
                   alpha_max :=alpha;
                end if;
              end findInterval;

            function solveOneNonlinearEquation
                "Solve f(u) = 0; f(u_min) and f(u_max) must have different signs"
                extends Modelica.Icons.Function;
                import Modelica.Utilities.Streams.error;

              input Real c1[:]
                  "[p] coefficients of denominator polynomials (c1[i]*p + 1)";
              input Real c2[:,2]
                  "[p^2, p] coefficients of denominator polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
              input Real u_min "Lower bound of search interval";
              input Real u_max "Upper bound of search interval";
              input Real tolerance=100*Modelica.Constants.eps
                  "Relative tolerance of solution u";
              output Real u "Value of independent variable so that f(u) = 0";

              protected
              constant Real eps=Modelica.Constants.eps "Machine epsilon";
              Real a=u_min "Current best minimum interval value";
              Real b=u_max "Current best maximum interval value";
              Real c "Intermediate point a <= c <= b";
              Real d;
              Real e "b - a";
              Real m;
              Real s;
              Real p;
              Real q;
              Real r;
              Real tol;
              Real fa "= f(a)";
              Real fb "= f(b)";
              Real fc;
              Boolean found=false;
            algorithm
              // Check that f(u_min) and f(u_max) have different sign
              fa := normalizationResidue(c1,c2,u_min);
              fb := normalizationResidue(c1,c2,u_max);
              fc := fb;
              if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
                error(
                  "The arguments u_min and u_max to solveOneNonlinearEquation(..)\n" +
                  "do not bracket the root of the single non-linear equation:\n" +
                  "  u_min  = " + String(u_min) + "\n" + "  u_max  = " + String(u_max)
                   + "\n" + "  fa = f(u_min) = " + String(fa) + "\n" +
                  "  fb = f(u_max) = " + String(fb) + "\n" +
                  "fa and fb must have opposite sign which is not the case");
              end if;

              // Initialize variables
              c := a;
              fc := fa;
              e := b - a;
              d := e;

              // Search loop
              while not found loop
                if abs(fc) < abs(fb) then
                  a := b;
                  b := c;
                  c := a;
                  fa := fb;
                  fb := fc;
                  fc := fa;
                end if;

                tol := 2*eps*abs(b) + tolerance;
                m := (c - b)/2;

                if abs(m) <= tol or fb == 0.0 then
                  // root found (interval is small enough)
                  found := true;
                  u := b;
                else
                  // Determine if a bisection is needed
                  if abs(e) < tol or abs(fa) <= abs(fb) then
                    e := m;
                    d := e;
                  else
                    s := fb/fa;
                    if a == c then
                      // linear interpolation
                      p := 2*m*s;
                      q := 1 - s;
                    else
                      // inverse quadratic interpolation
                      q := fa/fc;
                      r := fb/fc;
                      p := s*(2*m*q*(q - r) - (b - a)*(r - 1));
                      q := (q - 1)*(r - 1)*(s - 1);
                    end if;

                    if p > 0 then
                      q := -q;
                    else
                      p := -p;
                    end if;

                    s := e;
                    e := d;
                    if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*s*q) then
                      // interpolation successful
                      d := p/q;
                    else
                      // use bi-section
                      e := m;
                      d := e;
                    end if;
                  end if;

                  // Best guess value is defined as "a"
                  a := b;
                  fa := fb;
                  b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
                  fb := normalizationResidue(c1,c2,b);

                  if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
                    // initialize variables
                    c := a;
                    fc := fa;
                    e := b - a;
                    d := e;
                  end if;
                end if;
              end while;

              annotation (Documentation(info="<html>

<p>
This function determines the solution of <strong>one non-linear algebraic equation</strong> \"y=f(u)\"
in <strong>one unknown</strong> \"u\" in a reliable way. It is one of the best numerical
algorithms for this purpose. As input, the nonlinear function f(u)
has to be given, as well as an interval u_min, u_max that
contains the solution, i.e., \"f(u_min)\" and \"f(u_max)\" must
have a different sign. If possible, a smaller interval is computed by
inverse quadratic interpolation (interpolating with a quadratic polynomial
through the last 3 points and computing the zero). If this fails,
bisection is used, which always reduces the interval by a factor of 2.
The inverse quadratic interpolation method has superlinear convergence.
This is roughly the same convergence rate as a globally convergent Newton
method, but without the need to compute derivatives of the non-linear
function. The solver function is a direct mapping of the Algol 60 procedure
\"zero\" to Modelica, from:
</p>

<dl>
<dt> Brent R.P.:</dt>
<dd> <strong>Algorithms for Minimization without derivatives</strong>.
     Prentice Hall, 1973, pp. 58-59.</dd>
</dl>

</html>"));
            end solveOneNonlinearEquation;

            algorithm
               // Find interval for alpha
               (alpha_min, alpha_max) :=findInterval(c1, c2);

               // Compute alpha, so that abs(G(p)) = -3db
               alpha :=solveOneNonlinearEquation(
                c1,
                c2,
                alpha_min,
                alpha_max);
            end normalizationFactor;

            encapsulated function bandPassAlpha "Return alpha for band pass"
              extends Modelica.Icons.Function;

              import Modelica;
               input Real a "Coefficient of s^1";
               input Real b "Coefficient of s^0";
               input Modelica.Units.SI.AngularVelocity w
                "Bandwidth angular frequency";
               output Real alpha "Alpha factor to build up band pass";

            protected
              Real alpha_min;
              Real alpha_max;
              Real z_min;
              Real z_max;
              Real z;

              function residue "Residue of non-linear equation"
                extends Modelica.Icons.Function;
                input Real a;
                input Real b;
                input Real w;
                input Real z;
                output Real res;
              algorithm
                res := z^2 + (a*w*z/(1+z))^2 - (2+b*w^2)*z + 1;
              end residue;

            function solveOneNonlinearEquation
                "Solve f(u) = 0; f(u_min) and f(u_max) must have different signs"
                extends Modelica.Icons.Function;
                import Modelica.Utilities.Streams.error;

              input Real aa;
              input Real bb;
              input Real ww;
              input Real u_min "Lower bound of search interval";
              input Real u_max "Upper bound of search interval";
              input Real tolerance=100*Modelica.Constants.eps
                  "Relative tolerance of solution u";
              output Real u "Value of independent variable so that f(u) = 0";

              protected
              constant Real eps=Modelica.Constants.eps "Machine epsilon";
              Real a=u_min "Current best minimum interval value";
              Real b=u_max "Current best maximum interval value";
              Real c "Intermediate point a <= c <= b";
              Real d;
              Real e "b - a";
              Real m;
              Real s;
              Real p;
              Real q;
              Real r;
              Real tol;
              Real fa "= f(a)";
              Real fb "= f(b)";
              Real fc;
              Boolean found=false;
            algorithm
              // Check that f(u_min) and f(u_max) have different sign
              fa := residue(aa,bb,ww,u_min);
              fb := residue(aa,bb,ww,u_max);
              fc := fb;
              if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
                error(
                  "The arguments u_min and u_max to solveOneNonlinearEquation(..)\n" +
                  "do not bracket the root of the single non-linear equation:\n" +
                  "  u_min  = " + String(u_min) + "\n" + "  u_max  = " + String(u_max)
                   + "\n" + "  fa = f(u_min) = " + String(fa) + "\n" +
                  "  fb = f(u_max) = " + String(fb) + "\n" +
                  "fa and fb must have opposite sign which is not the case");
              end if;

              // Initialize variables
              c := a;
              fc := fa;
              e := b - a;
              d := e;

              // Search loop
              while not found loop
                if abs(fc) < abs(fb) then
                  a := b;
                  b := c;
                  c := a;
                  fa := fb;
                  fb := fc;
                  fc := fa;
                end if;

                tol := 2*eps*abs(b) + tolerance;
                m := (c - b)/2;

                if abs(m) <= tol or fb == 0.0 then
                  // root found (interval is small enough)
                  found := true;
                  u := b;
                else
                  // Determine if a bisection is needed
                  if abs(e) < tol or abs(fa) <= abs(fb) then
                    e := m;
                    d := e;
                  else
                    s := fb/fa;
                    if a == c then
                      // linear interpolation
                      p := 2*m*s;
                      q := 1 - s;
                    else
                      // inverse quadratic interpolation
                      q := fa/fc;
                      r := fb/fc;
                      p := s*(2*m*q*(q - r) - (b - a)*(r - 1));
                      q := (q - 1)*(r - 1)*(s - 1);
                    end if;

                    if p > 0 then
                      q := -q;
                    else
                      p := -p;
                    end if;

                    s := e;
                    e := d;
                    if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*s*q) then
                      // interpolation successful
                      d := p/q;
                    else
                      // use bi-section
                      e := m;
                      d := e;
                    end if;
                  end if;

                  // Best guess value is defined as "a"
                  a := b;
                  fa := fb;
                  b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
                  fb := residue(aa,bb,ww,b);

                  if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
                    // initialize variables
                    c := a;
                    fc := fa;
                    e := b - a;
                    d := e;
                  end if;
                end if;
              end while;

              annotation (Documentation(info="<html>

<p>
This function determines the solution of <strong>one non-linear algebraic equation</strong> \"y=f(u)\"
in <strong>one unknown</strong> \"u\" in a reliable way. It is one of the best numerical
algorithms for this purpose. As input, the nonlinear function f(u)
has to be given, as well as an interval u_min, u_max that
contains the solution, i.e., \"f(u_min)\" and \"f(u_max)\" must
have a different sign. If possible, a smaller interval is computed by
inverse quadratic interpolation (interpolating with a quadratic polynomial
through the last 3 points and computing the zero). If this fails,
bisection is used, which always reduces the interval by a factor of 2.
The inverse quadratic interpolation method has superlinear convergence.
This is roughly the same convergence rate as a globally convergent Newton
method, but without the need to compute derivatives of the non-linear
function. The solver function is a direct mapping of the Algol 60 procedure
\"zero\" to Modelica, from:
</p>

<dl>
<dt> Brent R.P.:</dt>
<dd> <strong>Algorithms for Minimization without derivatives</strong>.
     Prentice Hall, 1973, pp. 58-59.</dd>
</dl>

</html>"));
            end solveOneNonlinearEquation;

            algorithm
              assert( a^2/4 - b <= 0,  "Band pass transformation cannot be computed");
              z :=solveOneNonlinearEquation(a, b, w, 0, 1);
              alpha := sqrt(z);

              annotation (Documentation(info="<html>
<p>
A band pass with bandwidth \"w\" is determined from a low pass
</p>

<blockquote><pre>
1/(p^2 + a*p + b)
</pre></blockquote>

<p>
with the transformation
</p>

<blockquote><pre>
new(p) = (p + 1/p)/w
</pre></blockquote>

<p>
This results in the following derivation:
</p>

<blockquote><pre>
1/(p^2 + a*p + b) -> 1/( (p+1/p)^2/w^2 + a*(p + 1/p)/w + b )
                   = 1 /( ( p^2 + 1/p^2 + 2)/w^2 + (p + 1/p)*a/w + b )
                   = w^2*p^2 / (p^4 + 2*p^2 + 1 + (p^3 + p)a*w + b*w^2*p^2)
                   = w^2*p^2 / (p^4 + a*w*p^3 + (2+b*w^2)*p^2 + a*w*p + 1)
</pre></blockquote>

<p>
This 4th order transfer function shall be split in to two transfer functions of order 2 each
for numerical reasons. With the following formulation, the fourth order
polynomial can be represented (with the unknowns \"c\" and \"alpha\"):
</p>

<blockquote><pre>
g(p) = w^2*p^2 / ( (p*alpha)^2 + c*(p*alpha) + 1) * ( (p/alpha)^2 + c*(p/alpha) + 1)
     = w^2*p^2 / ( p^4 + c*(alpha + 1/alpha)*p^3 + (alpha^2 + 1/alpha^2 + c^2)*p^2
                                                 + c*(alpha + 1/alpha)*p + 1 )
</pre></blockquote>

<p>
Comparison of coefficients:
</p>

<blockquote><pre>
c*(alpha + 1/alpha) = a*w           -> c = a*w / (alpha + 1/alpha)
alpha^2 + 1/alpha^2 + c^2 = 2+b*w^2 -> equation to determine alpha

alpha^4 + 1 + a^2*w^2*alpha^4/(1+alpha^2)^2 = (2+b*w^2)*alpha^2
  or z = alpha^2
z^2 + a^2*w^2*z^2/(1+z)^2 - (2+b*w^2)*z + 1 = 0
</pre></blockquote>

<p>
Therefore the last equation has to be solved for \"z\" (basically, this means to compute
a real zero of a fourth order polynomial):
</p>

<blockquote><pre>
solve: 0 = f(z)  = z^2 + a^2*w^2*z^2/(1+z)^2 - (2+b*w^2)*z + 1  for \"z\"
           f(0)  = 1  &gt; 0
           f(1)  = 1 + a^2*w^2/4 - (2+b*w^2) + 1
                 = (a^2/4 - b)*w^2  &lt; 0
                 // since b - a^2/4 > 0 requirement for complex conjugate poles
-> 0 &lt; z &lt; 1
</pre></blockquote>

<p>
This function computes the solution of this equation and returns \"alpha = sqrt(z)\";
</p>

</html>"));
            end bandPassAlpha;
          end Utilities;
        end Filter;
      end Internal;
      annotation (
        Documentation(info="<html>
<p>
This package contains basic <strong>continuous</strong> input/output blocks
described by differential equations.
</p>

<p>
All blocks of this package can be initialized in different
ways controlled by parameter <strong>initType</strong>. The possible
values of initType are defined in
<a href=\"modelica://Modelica.Blocks.Types.Init\">Modelica.Blocks.Types.Init</a>:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Name</strong></td>
      <td><strong>Description</strong></td></tr>

  <tr><td><strong>Init.NoInit</strong></td>
      <td>no initialization (start values are used as guess values with fixed=false)</td></tr>

  <tr><td><strong>Init.SteadyState</strong></td>
      <td>steady state initialization (derivatives of states are zero)</td></tr>

  <tr><td><strong>Init.InitialState</strong></td>
      <td>Initialization with initial states</td></tr>

  <tr><td><strong>Init.InitialOutput</strong></td>
      <td>Initialization with initial outputs (and steady state of the states if possible)</td></tr>
</table>

<p>
For backward compatibility reasons the default of all blocks is
<strong>Init.NoInit</strong>, with the exception of Integrator and LimIntegrator
where the default is <strong>Init.InitialState</strong> (this was the initialization
defined in version 2.2 of the Modelica standard library).
</p>

<p>
In many cases, the most useful initial condition is
<strong>Init.SteadyState</strong> because initial transients are then no longer
present. The drawback is that in combination with a non-linear
plant, non-linear algebraic equations occur that might be
difficult to solve if appropriate guess values for the
iteration variables are not provided (i.e., start values with fixed=false).
However, it is often already useful to just initialize
the linear blocks from the Continuous blocks library in SteadyState.
This is uncritical, because only linear algebraic equations occur.
If Init.NoInit is set, then the start values for the states are
interpreted as <strong>guess</strong> values and are propagated to the
states with fixed=<strong>false</strong>.
</p>

<p>
Note, initialization with Init.SteadyState is usually difficult
for a block that contains an integrator
(Integrator, LimIntegrator, PI, PID, LimPID).
This is due to the basic equation of an integrator:
</p>

<blockquote><pre>
<strong>initial equation</strong>
   <strong>der</strong>(y) = 0;   // Init.SteadyState
<strong>equation</strong>
   <strong>der</strong>(y) = k*u;
</pre></blockquote>

<p>
The steady state equation leads to the condition that the input to the
integrator is zero. If the input u is already (directly or indirectly) defined
by another initial condition, then the initialization problem is <strong>singular</strong>
(has none or infinitely many solutions). This situation occurs often
for mechanical systems, where, e.g., u = desiredSpeed - measuredSpeed and
since speed is both a state and a derivative, it is always defined by
Init.InitialState or Init.SteadyState initialization.
</p>

<p>
In such a case, <strong>Init.NoInit</strong> has to be selected for the integrator
and an additional initial equation has to be added to the system
to which the integrator is connected. E.g., useful initial conditions
for a 1-dim. rotational inertia controlled by a PI controller are that
<strong>angle</strong>, <strong>speed</strong>, and <strong>acceleration</strong> of the inertia are zero.
</p>

</html>"),     Icon(graphics={Line(
              origin={0.061,4.184},
              points={{81.939,36.056},{65.362,36.056},{14.39,-26.199},{-29.966,
                  113.485},{-65.374,-61.217},{-78.061,-78.184}},
              color={95,95,95},
              smooth=Smooth.Bezier)}));
    end Continuous;

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

    package Math
    "Library of Real mathematical functions as input/output blocks"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.Package;

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

      block Product "Output product of the two inputs"
        extends Interfaces.SI2SO;

      equation
        y = u1*u2;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes the output <strong>y</strong>
as <em>product</em> of the two inputs <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
y = u1 * u2;
</pre></blockquote>

</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-100,60},{-40,60},{-30,40}}, color={0,0,127}),
              Line(points={{-100,-60},{-40,-60},{-30,-40}}, color={0,0,127}),
              Line(points={{50,0},{100,0}}, color={0,0,127}),
              Line(points={{-30,0},{30,0}}),
              Line(points={{-15,25.99},{15,-25.99}}),
              Line(points={{-15,-25.99},{15,25.99}}),
              Ellipse(lineColor={0,0,127}, extent={{-50,-50},{50,50}})}));
      end Product;

      block Division "Output first input divided by second input"
        extends Interfaces.SI2SO;

      equation
        y = u1/u2;
        annotation (
          Documentation(info="<html>
<p>
This block computes the output <strong>y</strong>
by <em>dividing</em> the two inputs <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
y = u1 / u2;
</pre></blockquote>

</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-100,60},{-60,60},{0,0}}, color={0,0,127}),
              Line(points={{-100,-60},{-60,-60},{0,0}}, color={0,0,127}),
              Ellipse(lineColor={0,0,127}, extent={{-50,-50},{50,50}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(points={{50,0},{100,0}}, color={0,0,127}),
              Line(points={{-30,0},{30,0}}),
              Ellipse(fillPattern=FillPattern.Solid, extent={{-5,20},{5,30}}),
              Ellipse(fillPattern=FillPattern.Solid, extent={{-5,-30},{5,-20}}),
              Text(
                extent={{-60,90},{90,50}},
                textColor={128,128,128},
                textString="u1 / u2")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={         Line(points={{50,0},{100,0}},
                color={0,0,255})}));
      end Division;

      block Min "Pass through the smallest signal"
        extends Interfaces.SI2SO;
      equation
        y = min(u1, u2);
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={Text(
                extent={{-90,36},{90,-36}},
                textColor={160,160,164},
                textString="min()")}), Documentation(info="<html>
<p>
This block computes the output <strong>y</strong> as <em>minimum</em> of
the two Real inputs <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
y = <strong>min</strong> ( u1 , u2 );
</pre></blockquote>
</html>"));
      end Min;
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

    package Types
    "Library of constants, external objects and types with choices, especially to build menus"
      extends Modelica.Icons.TypesPackage;

        type Init = enumeration(
          NoInit
            "No initialization (start values are used as guess values with fixed=false)",
          SteadyState
            "Steady state initialization (derivatives of states are zero)",
          InitialState "Initialization with initial states",
          InitialOutput
            "Initialization with initial outputs (and steady state of the states if possible)")
        "Enumeration defining initialization of a block" annotation (Evaluate=true,
        Documentation(info="<html>
  <p>The following initialization alternatives are available:</p>
  <dl>
    <dt><code><strong>NoInit</strong></code></dt>
      <dd>No initialization (start values are used as guess values with <code>fixed=false</code>)</dd>
    <dt><code><strong>SteadyState</strong></code></dt>
      <dd>Steady state initialization (derivatives of states are zero)</dd>
    <dt><code><strong>InitialState</strong></code></dt>
      <dd>Initialization with initial states</dd>
    <dt><code><strong>InitialOutput</strong></code></dt>
      <dd>Initialization with initial outputs (and steady state of the states if possible)</dd>
  </dl>
</html>"));

      type AnalogFilter = enumeration(
          CriticalDamping "Filter with critical damping",
          Bessel "Bessel filter",
          Butterworth "Butterworth filter",
          ChebyshevI "Chebyshev I filter")
        "Enumeration defining the method of filtering" annotation (Evaluate=true);

      type FilterType = enumeration(
          LowPass "Low pass filter",
          HighPass "High pass filter",
          BandPass "Band pass filter",
          BandStop "Band stop / notch filter")
        "Enumeration of analog filter types (low, high, band pass or band stop filter)"
        annotation (Evaluate=true);
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

  package Icons "Icons for Math"
    extends Modelica.Icons.IconsPackage;

    partial function AxisLeft
      "Basic icon for mathematical function with y-axis on left side"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,-80},{-80,68}}, color={192,192,192}),
            Polygon(
              points={{-80,90},{-88,68},{-72,68},{-80,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              textColor={0,0,255})}),
        Documentation(info="<html>
<p>
Icon for a mathematical function, consisting of an y-axis on the left side.
It is expected, that an x-axis is added and a plot of the function.
</p>
</html>"));
    end AxisLeft;

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

  function cos "Cosine"
    extends Modelica.Math.Icons.AxisLeft;
    input Modelica.Units.SI.Angle u "Independent variable";
    output Real y "Dependent variable y=cos(u)";

  external "builtin" y = cos(u);
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
          Line(points={{-80,80},{-74.4,78.1},{-68.7,72.3},{-63.1,63},{-56.7,48.7},
                {-48.6,26.6},{-29.3,-32.5},{-22.1,-51.7},{-15.7,-65.3},{-10.1,-73.8},
                {-4.42,-78.8},{1.21,-79.9},{6.83,-77.1},{12.5,-70.6},{18.1,-60.6},
                {24.5,-45.7},{32.6,-23},{50.3,31.3},{57.5,50.7},{63.9,64.6},{69.5,
                73.4},{75.2,78.6},{80,80}}),
          Text(
            extent={{-36,82},{36,34}},
            textColor={192,192,192},
            textString="cos")}),
      Documentation(info="<html>
<p>
This function returns y = cos(u), with -&infin; &lt; u &lt; &infin;:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/cos.png\">
</p>
</html>"));
  end cos;

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

  function asinh "Inverse of sinh (area hyperbolic sine)"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u "Independent variable";
    output Real y "Dependent variable y=asinh(u)";

  algorithm
    y := Modelica.Math.log(u + sqrt(u*u + 1));
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
          Line(points={{-80,-80},{-56.7,-68.4},{-39.8,-56.8},{-26.9,-44.7},{-17.3,
                -32.4},{-9.25,-19},{9.25,19},{17.3,32.4},{26.9,44.7},{39.8,56.8},
                {56.7,68.4},{80,80}}),
          Text(
            extent={{-90,80},{-6,26}},
            textColor={192,192,192},
            textString="asinh")}),
      Documentation(info="<html>
<p>
The function returns the area hyperbolic sine of its
input argument u. This inverse of sinh(..) is unique
and there is no restriction on the input argument u of
asinh(u) (-&infin; &lt; u &lt; &infin;):
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/asinh.png\">
</p>
</html>"));
  end asinh;

  function log "Natural (base e) logarithm (u shall be > 0)"
    extends Modelica.Math.Icons.AxisLeft;
    input Real u "Independent variable";
    output Real y "Dependent variable y=ln(u)";

  external "builtin" y = log(u);
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
          Line(points={{-80,-80},{-79.2,-50.6},{-78.4,-37},{-77.6,-28},{-76.8,-21.3},
                {-75.2,-11.4},{-72.8,-1.31},{-69.5,8.08},{-64.7,17.9},{-57.5,28},
                {-47,38.1},{-31.8,48.1},{-10.1,58},{22.1,68},{68.7,78.1},{80,80}}),
          Text(
            extent={{-6,-24},{66,-72}},
            textColor={192,192,192},
            textString="log")}),
      Documentation(info="<html>
<p>
This function returns y = log(10) (the natural logarithm of u),
with u &gt; 0:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/log.png\">
</p>
</html>"));
  end log;
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

    package Streams "Read from files and write to files"
      extends Modelica.Icons.FunctionsPackage;

      pure function error "Print error message and cancel all actions - in case of an unrecoverable error"
        extends Modelica.Icons.Function;
        input String string "String to be printed to error message window";
        external "C" ModelicaError(string) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaUtilities.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Streams.<strong>error</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
In case of an unrecoverable error (i.e., if the solver is unable to recover from the error),
print the string \"string\" as error message and cancel all actions.
This function is semantically equivalent with the built-in function <strong>assert</strong> if called with the (default) <strong>AssertionLevel.error</strong>.
Line breaks are characterized by \"\\n\" in the string.
</p>
<h4>Example</h4>
<blockquote><pre>
Streams.error(\"x (= \" + String(x) + \")\\nhas to be in the range 0 .. 1\");
</pre></blockquote>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams\">Streams</a>,
<a href=\"modelica://Modelica.Utilities.Streams.print\">Streams.print</a>,
<a href=\"modelica://ModelicaReference.Operators.'assert()'\">ModelicaReference.Operators.'assert()'</a>
<a href=\"modelica://ModelicaReference.Operators.'String()'\">ModelicaReference.Operators.'String()'</a>
</p>
</html>"));
      end error;
      annotation (
        Documentation(info="<html>
<h4>Library content</h4>
<p>
Package <strong>Streams</strong> contains functions to input and output strings
to a message window or on files, as well as reading matrices from file
and writing matrices to file. Note that a string is interpreted
and displayed as html text (e.g., with print(..) or error(..))
if it is enclosed with the Modelica html quotation, e.g.,
</p>
<blockquote><p>
string = \"&lt;html&gt; first line &lt;br&gt; second line &lt;/html&gt;\".
</p></blockquote>
<p>
It is a quality of implementation, whether (a) all tags of html are supported
or only a subset, (b) how html tags are interpreted if the output device
does not allow to display formatted text.
</p>
<p>
In the table below an example call to every function is given:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><th><strong><em>Function/type</em></strong></th><th><strong><em>Description</em></strong></th></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.print\">print</a>(string)<br>
          <a href=\"modelica://Modelica.Utilities.Streams.print\">print</a>(string,fileName)</td>
      <td> Print string \"string\" or vector of strings to message window or on
           file \"fileName\".</td>
  </tr>
  <tr><td>stringVector =
         <a href=\"modelica://Modelica.Utilities.Streams.readFile\">readFile</a>(fileName)</td>
      <td> Read complete text file and return it as a vector of strings.</td>
  </tr>
  <tr><td>(string, endOfFile) =
         <a href=\"modelica://Modelica.Utilities.Streams.readLine\">readLine</a>(fileName, lineNumber)</td>
      <td>Returns from the file the content of line lineNumber.</td>
  </tr>
  <tr><td>lines =
         <a href=\"modelica://Modelica.Utilities.Streams.countLines\">countLines</a>(fileName)</td>
      <td>Returns the number of lines in a file.</td>
  </tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.error\">error</a>(string)</td>
      <td> Print error message \"string\" to message window
           and cancel all actions</td>
  </tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.close\">close</a>(fileName)</td>
      <td> Close file if it is still open. Ignore call if
           file is already closed or does not exist. </td>
  </tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">readMatrixSize</a>(fileName, matrixName)</td>
      <td> Read dimensions of a Real matrix from a MATLAB MAT file. </td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">readRealMatrix</a>(fileName, matrixName, nrow, ncol)</td>
      <td> Read a Real matrix from a MATLAB MAT file. </td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.writeRealMatrix\">writeRealMatrix</a>(fileName, matrixName, matrix, append, format)</td>
      <td> Write Real matrix to a MATLAB MAT file. </td></tr>
</table>
<p>
Use functions <strong>scanXXX</strong> from package
<a href=\"modelica://Modelica.Utilities.Strings\">Strings</a>
to parse a string.
</p>
<p>
If Real, Integer or Boolean values shall be printed
or used in an error message, they have to be first converted
to strings with the builtin operator
<a href=\"modelica://ModelicaReference.Operators.'String()'\">ModelicaReference.Operators.'String()'</a>(...).
Example:
</p>
<blockquote><pre>
<strong>if</strong> x &lt; 0 <strong>or</strong> x &gt; 1 <strong>then</strong>
   Streams.error(\"x (= \" + String(x) + \") has to be in the range 0 .. 1\");
<strong>end if</strong>;
</pre></blockquote>
</html>"));
    end Streams;

    package Strings "Operations on strings"
      extends Modelica.Icons.FunctionsPackage;

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

      type AngularVelocity = Real (
          final quantity="AngularVelocity",
          final unit="rad/s");

      type Velocity = Real (final quantity="Velocity", final unit="m/s");

      type Acceleration = Real (final quantity="Acceleration", final unit="m/s2");

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

      type AmountOfSubstance = Real (
          final quantity="AmountOfSubstance",
          final unit="mol",
          min=0);

      type Molality = Real (final quantity="Molality", final unit="mol/kg");

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

      package Air
        import Physiolibrary.Media.Substances.*;
        extends Interfaces.PartialMedium(
           ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
           reducedX = false,
           singleState = false,
           substanceNames={"O2","CO2","H2O","N2"},
           reference_X=cat(1, Conc .* C2X, {1 - (Conc * C2X)}),
           SpecificEnthalpy(start=0, nominal=1e3),
           Density(start=1.0, nominal=1.0),
           AbsolutePressure(start=1.0e5, nominal=1.0e5),
           Temperature(min=273.15, max=320.15, start=298.15, nominal=298.15),
           MassFlowRate(nominal=1e-3));

  protected
        package stateOfMatter = Chemical.Interfaces.IdealGas
          "Substances model to translate data into substance properties";

        constant stateOfMatter.SubstanceData substanceData[nS]={
            Chemical.Substances.Oxygen_gas(),
            Chemical.Substances.CarbonDioxide_gas(),
            Chemical.Substances.Water_gas(),
            Chemical.Substances.Nitrogen_gas()} "Definition of the substances";

        constant Modelica.Units.SI.MoleFraction Conc[nS-1]={0.21,0.0004,0.02}
          "sum(*) = 1";

        constant Real C2X[nS-1] = aMM[1:nS-1] ./ (Conc * aMM[1:nS-1] + (1 - sum(Conc)) * aMM[nS]) "Conc to mass fraction coefficient";

        constant Real aMM[nS] = ones(nS) ./ stateOfMatter.specificAmountOfParticles(substanceData, T=298.15, p=101325) "Average molar mass of substance particle";

  public
        redeclare connector extends SubstancesPort
         Chemical.Interfaces.SubstancePort_a O2 "Gaseous oxygen molecule";
         Chemical.Interfaces.SubstancePort_a CO2 "Gaseous hydrogen molecule";
         Chemical.Interfaces.SubstancePort_a H2O "Gaseous H2O molecule";
         Chemical.Interfaces.SubstancePort_a N2 "Gaseaous nitrogen molecule";
        end SubstancesPort;

        redeclare replaceable model extends SubstancesDecomposition "Just because Modelica in today version cannot work properly with nested connectors"

        Chemical.Interfaces.SubstancePort_a O2 annotation (Placement(transformation(extent={{90,50},{110,70}})));
        Chemical.Interfaces.SubstancePort_a CO2 annotation (Placement(transformation(extent={{90,90},{110,110}})));
        Chemical.Interfaces.SubstancePort_a N2 annotation (Placement(transformation(extent={{90,-70},{110,-50}}), iconTransformation(extent={{90,-70},{110,-50}})));
        Chemical.Interfaces.SubstancePort_a H2O annotation (Placement(transformation(extent={{90,-110},{110,-90}})));
        equation
        connect(O2, substances.O2) annotation (Line(points={{100,60},{-72,60},{-72,0},{-100,0}},      color={158,66,200}));
        connect(CO2, substances.CO2) annotation (Line(points={{100,100},{22,100},{22,80},{-76,80},{-76,0},{-100,0},{-100,0}},     color={158,66,200}));
        connect(N2, substances.N2) annotation (Line(points={{100,-60},{88,-60},{88,-46},{-76,-46},{-76,0},{-100,0}},   color={158,66,200}));
        connect(H2O, substances.H2O) annotation (Line(points={{100,-100},{-82,-100},{-82,0},{-100,0}},      color={158,66,200}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
        end SubstancesDecomposition;


        redeclare replaceable model extends ChemicalSolution
    protected
          Real I = 0 "mole-fraction-based ionic strength";
           Modelica.Units.SI.Molality NpM[nS] "Amount of substance particles per mass of substance";
           Modelica.Units.SI.MoleFraction x_baseMolecule[nS] "Mole fraction of free base molecule of substance";


        //initial equation
        //  substanceMasses = startSubstanceMasses;
        equation
          v=0 "electric potential is not used without external flows of charge";

          NpM = stateOfMatter.specificAmountOfParticles(substanceData,T=T,p=p);

          x_baseMolecule = X.*stateOfMatter.specificAmountOfFreeBaseMolecule(substanceData,T=T,p=p)./(X*NpM);

          T = stateOfMatter.solution_temperature(
              substanceData,
              h,
              X,
              p);


          substances.O2.u =stateOfMatter.electroChemicalPotentialPure(
              Substances.O2_g,
              T,
              p,
              v,
              I) + Modelica.Constants.R*T*log(x_baseMolecule[
            i("O2")]);
          substances.CO2.u =stateOfMatter.electroChemicalPotentialPure(
              Substances.CO2_g,
              T,
              p,
              v,
              I) + Modelica.Constants.R*T*log(x_baseMolecule[
            i("CO2")]);
          substances.H2O.u =stateOfMatter.electroChemicalPotentialPure(
              Substances.H2O_g,
              T,
              p,
              v,
              I) + Modelica.Constants.R*T*log(x_baseMolecule[
            i("H2O")]);
          substances.N2.u =stateOfMatter.electroChemicalPotentialPure(
              Substances.N2_g,
              T,
              p,
              v,
              I) + Modelica.Constants.R*T*log(x_baseMolecule[
            i("N2")]);
          substances.O2.h_outflow = stateOfMatter.molarEnthalpy( Substances.O2_g, T, p, v, I);

          substances.CO2.h_outflow = stateOfMatter.molarEnthalpy( Substances.CO2_g, T, p, v, I);
          substances.H2O.h_outflow = stateOfMatter.molarEnthalpy( Substances.H2O_g, T, p, v, I);
          substances.N2.h_outflow = stateOfMatter.molarEnthalpy( Substances.N2_g, T, p, v, I);

          enthalpyFromSubstances =
           substances.O2.q * actualStream(substances.O2.h_outflow) +
           substances.CO2.q * actualStream(substances.CO2.h_outflow) +
           substances.H2O.q * actualStream(substances.H2O.h_outflow) +
           substances.N2.q * actualStream(substances.N2.h_outflow)
            "enthalpy from substances";


          massFlows[i("O2")] = substances.O2.q*Substances.O2_g.MolarWeight;
          massFlows[i("CO2")] = substances.CO2.q*Substances.CO2_g.MolarWeight;
          massFlows[i("H2O")] = substances.H2O.q*Substances.H2O_g.MolarWeight;
          massFlows[i("N2")] = substances.N2.q*Substances.N2_g.MolarWeight;

        end ChemicalSolution;

        redeclare function extends specificEnthalpies_TpvI
        algorithm
             specificEnthalpy:=stateOfMatter.specificEnthalpy(
                substanceData,
                T,p,v,I);
        end specificEnthalpies_TpvI;

  public
        redeclare replaceable record extends ThermodynamicState
          "A selection of variables that uniquely defines the thermodynamic state"
          extends Modelica.Icons.Record;
          AbsolutePressure p "Absolute pressure of medium";
          Temperature T "Temperature of medium";
          MassFraction X[nS] "Mass fractions of substances";
          annotation (Documentation(info="<html>

</html>"));
        end ThermodynamicState;

        replaceable function electrochemicalPotentials_pTXvI
          input Modelica.Units.SI.Pressure p;
          input Modelica.Units.SI.Temperature T;
          input Modelica.Units.SI.MoleFraction x_baseMolecule[nS] "Free mole fraction of substance base molecule";
          input Modelica.Units.SI.ElectricPotential electricPotential=0;
          input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0;
          output Modelica.Units.SI.ChemicalPotential u[nS];
    protected
          Real a[nS];
          Modelica.Units.SI.ChargeNumberOfIon z[nS];
        algorithm
          a := stateOfMatter.activityCoefficient(
              substanceData,
              T,
              p,
              electricPotential,
              moleFractionBasedIonicStrength) .* x_baseMolecule;
          z := stateOfMatter.chargeNumberOfIon(
              substanceData,
              T,
              p,
              electricPotential,
              moleFractionBasedIonicStrength);
          u := stateOfMatter.chemicalPotentialPure(
              substanceData,
              T,
              p,
              electricPotential,
              moleFractionBasedIonicStrength) .+ Modelica.Constants.R*T*log(a) .+ z*Modelica.Constants.F*
            electricPotential;
        end electrochemicalPotentials_pTXvI;

        redeclare replaceable function extends setState_pTX
        algorithm
          state.T := T;
          state.p := p;
          state.X := X;
        end setState_pTX;

        redeclare replaceable function extends setState_phX
          "Return thermodynamic state as function of p, h and composition X or Xi"
        algorithm
          state.p := p;
          state.X := X;
          state.T := stateOfMatter.solution_temperature(
              substanceData,
              h,
              X,
              p);
        end setState_phX;

         redeclare replaceable function extends density
         algorithm
          d := 1/(state.X*stateOfMatter.specificVolume(
              substanceData,
              state.T,
              state.p));
         end density;

        redeclare replaceable function extends specificEnthalpy
        algorithm
          h := state.X * stateOfMatter.specificEnthalpy(
              substanceData,
              state.T,
              state.p);
          /*, electricPotential, moleFractionBasedIonicStrength*/
        end specificEnthalpy;

        redeclare replaceable function extends temperature
        algorithm
          T := state.T;
        end temperature;

        redeclare replaceable function extends pressure
        algorithm
          p := state.p;
        end pressure;

        function X "To set mass fractions"
          input Types.AmountOfSubstance
              tO2 = 0.21,
              tCO2 = 0.0003,
              tH2O = 0.06,
              tN2 = 1-tO2-tCO2-tH2O;
          output Types.MassFraction X[nS];
    protected
          Types.Mass tm;
        algorithm
          tm :=tO2*O2_g.MolarWeight + tCO2*CO2_g.MolarWeight + tH2O*H2O_g.MolarWeight +
            tN2*N2_g.MolarWeight;
          X[i("O2")] := (tO2*O2_g.MolarWeight)/tm;
          X[i("CO2")] := (tCO2*CO2_g.MolarWeight)/tm;
          X[i("H2O")] := (tH2O*H2O_g.MolarWeight)/tm;
          X[i("N2")] := (tN2*N2_g.MolarWeight)/tm;
        end X;



        annotation (Documentation(revisions="<html>
<p><i>2021</i></p>
<p>Marek Matejak, http://www.physiolib.com </p>
<p>All rights reserved. </p>
</html>"));
      end Air;

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

      constant Chemical.Interfaces.IdealGas.SubstanceData O2_g=
            Chemical.Substances.Oxygen_gas();

        constant Chemical.Interfaces.IdealGas.SubstanceData CO2_g=
            Chemical.Substances.CarbonDioxide_gas();

        constant Chemical.Interfaces.IdealGas.SubstanceData H2O_g=
            Chemical.Substances.Water_gas();

        constant Chemical.Interfaces.IdealGas.SubstanceData N2_g=
            Chemical.Substances.Nitrogen_gas();

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

      model VolumePump "Prescribed mass flow"
        extends Icons.Pump;
        extends Physiolibrary.Fluid.Interfaces.OnePort;
        extends Physiolibrary.Fluid.Interfaces.ConditionalVolumeFlow;
      equation
        volumeFlowRate = q;
        annotation (
          Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-150, -90}, {150, -50}}, lineColor = {127, 0, 0}, textString = "%name")}),
          Documentation(revisions = "<html>
<table>
<tr>
<td>Author:</td>
<td>Marek Matejak</td>
</tr>
<td>Web:</td>
<td>http://www.physiolib.com</td>
</tr>
<tr>
<td>Date of:</td>
<td>october 2017-2018</td>
</tr>
</table>
</html>", info = "<html>
<p><font style=\"font-size: 9pt; \">This element needs to be connected only to next hydraulic elements, which contain calculation of hydraulic pressure in connector. It is because equation contains only </font><b><font style=\"font-size: 9pt; \">hydraulic volume flow</font></b><font style=\"font-size: 9pt; \"> variable, which is set to value of input signal variable. </font></p>
</html>"));
      end VolumePump;
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

      partial model ConditionalVolumeFlow "Input of solution volume flow vs. parametric solution volume flow"
        parameter Boolean useSolutionFlowInput = false "=true, if solution flow input is used instead of parameter SolutionFlow" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));
        parameter Physiolibrary.Types.VolumeFlowRate SolutionFlow = 0 "Mass flow of solution if useSolutionFlowInput=false" annotation (
          HideResult = not useSolutionFlowInput,
          Dialog(enable = not useSolutionFlowInput));
        Physiolibrary.Types.RealIO.VolumeFlowRateInput solutionFlow(start = SolutionFlow) = q if useSolutionFlowInput annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 40}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, 70})));
        Physiolibrary.Types.VolumeFlowRate q "Current solution flow";
      equation
        if not useSolutionFlowInput then
          q = SolutionFlow;
        end if;
      end ConditionalVolumeFlow;

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

      model FlowMeasure "Volumetric flow between ports"
        extends Physiolibrary.Fluid.Interfaces.OnePort;
        // extends Icons.FlowMeasure;
        extends Modelica.Icons.RoundSensor;
        Physiolibrary.Types.RealIO.MassFlowRateOutput massFlow
        "Actual mass flow rate"                                                        annotation (
          Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 270, origin = {0, -60}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {0, 120})));
        Physiolibrary.Types.VolumeFlowRate volumeInflowRate;
        Physiolibrary.Types.VolumeFlowRate volumeOutflowRate;
      protected
        Medium.ThermodynamicState state_inflow "state for medium inflowing through q_in";
        Medium.ThermodynamicState state_outflow "state for medium outflowing through q_out";
        Modelica.Units.SI.Density density_inflow;
        Modelica.Units.SI.Density density_outflow;
      equation
        q_out.p = q_in.p;
        massFlow = q_in.m_flow;
      // medium states
        state_inflow = Medium.setState_phX(q_in.p, inStream(q_in.h_outflow), inStream(q_in.Xi_outflow));
        state_outflow = Medium.setState_phX(q_out.p, inStream(q_out.h_outflow), inStream(q_out.Xi_outflow));
        density_inflow = Medium.density(state_inflow);
        density_outflow = Medium.density(state_outflow);
        volumeInflowRate = massFlow / density_inflow;
        volumeOutflowRate = massFlow / density_outflow;
        annotation (
          Documentation(revisions = "<html>
	<p><i>2009-2018</i></p>
	<p>Marek Matejak, marek@matfyz.cz </p>
	</html>"),
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-25, -11}, {34, -70}}, lineColor = {0, 0, 0}, textString = "V'")}));
      end FlowMeasure;

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
    end Sensors;

    package Sources
      extends Modelica.Icons.SourcesPackage;

      model PressureSource "Prescribed pressure at port with unlimited mass"
        parameter Modelica.Units.SI.Pressure pressure_start = system.p_ambient "Initial pressure" annotation (
          Dialog(enable = not usePressureInput, group = "Initialization"));
        extends Physiolibrary.Fluid.Interfaces.CompositionSetup;

        parameter Boolean usePressureInput = false "=true, if pressure input is used" annotation (
          Evaluate = true,
          HideResult = true,
          choices(checkBox = true),
          Dialog(group = "Conditional inputs"));

        Physiolibrary.Types.RealIO.PressureInput pressure(start = pressure_start) = p if usePressureInput "Pressure" annotation (
          Placement(transformation(extent = {{-120, -20}, {-80, 20}})));
        Physiolibrary.Fluid.Interfaces.FluidPort_a y(redeclare package Medium
          =                                                                     Medium) "PressureFlow output connectors" annotation (
          Placement(transformation(extent = {{84, -16}, {116, 16}})));

      protected
        Physiolibrary.Types.Pressure p;
        Modelica.Units.SI.SpecificEnthalpy h = Medium.specificEnthalpy_pTX(p, temperature_start, x_mass_start) "Fluid enthalphy";
      equation
        if not usePressureInput then
          p = pressure_start;
        end if;
        y.p = p;
        y.h_outflow = h;
        y.Xi_outflow = x_mass_start[1:Medium.nXi];
        y.C_outflow = C_start;
        annotation (
          Documentation(info = "<html>
        <p>Model has a vector of continuous Real input signals as pressures for vector of pressure-flow connectors. </p>
        <p>Usage in tests: Set defaul volume&gt;0 and try to set STEADY in instances to &quot;false&quot;!</p>
        </html>", revisions = "<html>
        <p><i>2009-2018</i></p>
        <p>Marek Matejak, marek@matfyz.cz </p>
        </html>"),
          Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, pattern = LinePattern.None, fillColor = {170, 255, 255}, fillPattern = FillPattern.Backward), Text(extent = {{0, 0}, {-100, -100}}, lineColor = {0, 0, 0}, textString = "P"), Line(points = {{-100, 0}, {56, 0}}, color = {191, 0, 0}, thickness = 0.5), Polygon(points = {{38, -20}, {38, 20}, {78, 0}, {38, -20}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 150}, {150, 110}}, textString = "%name", lineColor = {0, 0, 255})}));
      end PressureSource;
    end Sources;
    annotation (
      Documentation(info = "<html>
<p>The main usage of this fluid domain is modeling of the cardio-vascular, respiratory and lymhpatic system in human physiology. And because there are no extreme thermodynamic conditions, the system can be really simple &mdash;it is only necessary to model conditions for ideal gases, for incompressible liquids, at normal liquid temperatures and with relative pressure 5-20kPa. This boring thermodynamic state leads to the very simple blocks of resistance,  pressure, volumetric flow, inertia and finally the block of blood accumulation in elastic comparments.</p>
</html>"));
  end Fluid;

  package Icons "Icons for physiological models"
    extends Modelica.Icons.Package;

    class RespiratoryCenter
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/respiracniCentrum.png")}));
    end RespiratoryCenter;

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

    class IdealValve
      annotation (
        Icon(graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/ideal_valve.svg")}));
    end IdealValve;

    class Pump
      annotation (
        Icon(graphics={  Bitmap(extent = {{-100, -100}, {100, 100}}, fileName = "modelica://Physiolibrary/Resources/Icons/pump.svg")}));
    end Pump;
    annotation (
      Documentation(revisions = ""));
  end Icons;

  package Types "Physiological units with nominals"
    extends Modelica.Icons.Package;

    package Constants
      extends Modelica.Icons.SourcesPackage;

      block FrequencyConst "Constant signal of type Frequency"
        parameter Types.Frequency k "Constant Frequency output value";
        RealIO.FrequencyOutput y "Frequency constant" annotation (
          Placement(transformation(extent = {{40, -10}, {60, 10}}), iconTransformation(extent = {{40, -10}, {60, 10}})));
      equation
        y = k;
        annotation (
          defaultComponentName = "frequency",
          Diagram(coordinateSystem(extent = {{-40, -40}, {40, 40}})),
          Icon(coordinateSystem(extent = {{-40, -40}, {40, 40}}, preserveAspectRatio = false), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, radius = 10, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -44}, {100, -64}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "%name"), Text(extent = {{-40, 10}, {40, -10}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "Const")}));
      end FrequencyConst;

      block HydraulicConductanceConst "Constant signal of type HydraulicConductance"
        parameter Types.HydraulicConductance k "Constant HydraulicConductance output value";
        RealIO.HydraulicConductanceOutput y "HydraulicConductance constant" annotation (
          Placement(transformation(extent = {{40, -10}, {60, 10}}), iconTransformation(extent = {{40, -10}, {60, 10}})));
      equation
        y = k;
        annotation (
          defaultComponentName = "hydraulicConductance",
          Diagram(coordinateSystem(extent = {{-40, -40}, {40, 40}})),
          Icon(coordinateSystem(extent = {{-40, -40}, {40, 40}}, preserveAspectRatio = false), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, radius = 10, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -44}, {100, -64}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "%name"), Text(extent = {{-40, 10}, {40, -10}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "Const")}));
      end HydraulicConductanceConst;

      block VolumeConst "Constant signal of type Volume"
        parameter Types.Volume k "Constant Volume output value";
        RealIO.VolumeOutput y "Volume constant" annotation (
          Placement(transformation(extent = {{40, -10}, {60, 10}}), iconTransformation(extent = {{40, -10}, {60, 10}})));
      equation
        y = k;
        annotation (
          defaultComponentName = "volume",
          Diagram(coordinateSystem(extent = {{-40, -40}, {40, 40}})),
          Icon(coordinateSystem(extent = {{-40, -40}, {40, 40}}, preserveAspectRatio = false), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, radius = 10, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -44}, {100, -64}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "%name"), Text(extent = {{-40, 10}, {40, -10}}, lineColor = {0, 0, 0}, fillColor = {236, 236, 236}, fillPattern = FillPattern.Solid, textString = "Const")}));
      end VolumeConst;
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

      connector VolumeFlowRateInput = input VolumeFlowRate "input VolumeFlowRate as connector" annotation (
        defaultComponentName = "volumeflowrate",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type VolumeFlowRate.
            </p>
            </html>"));

      connector VolumeFlowRateOutput = output VolumeFlowRate "output VolumeFlowRate as connector" annotation (
        defaultComponentName = "volumeflowrate",
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

      connector ElectricPotentialOutput = output ElectricPotential "output ElectricPotential as connector" annotation (
        defaultComponentName = "electricpotential",
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

      connector FrequencyOutput = output Frequency "output Frequency as connector" annotation (
        defaultComponentName = "frequency",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Frequency.
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

      connector HydraulicConductanceOutput = output HydraulicConductance "output HydraulicConductance as connector" annotation (
        defaultComponentName = "hydraulicconductance",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{30, 110}, {30, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
          <p>
          Connector with one output signal of type Real.
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

      connector SpecificEnthalpyInput = input SpecificEnthalpy "input SpecificEnthalpy as connector" annotation (
        defaultComponentName = "specificEnthalpy",
        Icon(graphics={  Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2)),
        Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-10, 85}, {-10, 60}}, lineColor = {0, 0, 127}, textString = "%name")}),
        Documentation(info = "<html>
            <p>
            Connector with one input signal of type SpecificEnthalpy.
            </p>
            </html>"));
    end RealIO;

    type Frequency = Modelica.Units.SI.Frequency(displayUnit = "1/min");

    type Mass = Modelica.Units.SI.Mass(displayUnit = "g", nominal = 1e-3, min = 0, max = Modelica.Constants.inf);

    type MassFraction = Modelica.Units.SI.MassFraction(nominal = 0.1, min = ModelicaServices.Machine.small, max = Modelica.Constants.inf);

    type MassFlowRate = Modelica.Units.SI.MassFlowRate(displayUnit = "mg/min", nominal = 0.001);

    type Density = Modelica.Units.SI.Density(displayUnit = "kg/l", nominal = 1e-3);

    type Pressure = Modelica.Units.SI.Pressure(displayUnit = "mmHg", nominal = 1e5);

    type Volume = Modelica.Units.SI.Volume(displayUnit = "ml", nominal = 1e-6, min = 0, max = Modelica.Constants.inf);

    type VolumeFlowRate = Modelica.Units.SI.VolumeFlowRate(displayUnit = "ml/min", nominal = 1e-6 / 60);

    replaceable type Concentration = Modelica.Units.SI.Concentration(displayUnit = "mmol/l", min = ModelicaServices.Machine.small, max = Modelica.Constants.inf) constrainedby Real;

    replaceable type AmountOfSubstance = Modelica.Units.SI.AmountOfSubstance(displayUnit = "mmol", min = 0, max = Modelica.Constants.inf) constrainedby Real;

    type MassConcentration = Modelica.Units.SI.MassConcentration(displayUnit = "mg/l", nominal = 1e-3, min = ModelicaServices.Machine.small, max = Modelica.Constants.inf);

    type Temperature = Modelica.Units.SI.Temperature(displayUnit = "degC", nominal = 1, min = 0);

    type HeatFlowRate = Modelica.Units.SI.HeatFlowRate(displayUnit = "kcal/min", nominal = 4186.8 / 60);

    type SpecificEnthalpy = Modelica.Units.SI.SpecificEnthalpy(displayUnit = "kcal/kg", nominal = 1e5);

    type ElectricPotential = Modelica.Units.SI.ElectricPotential(displayUnit = "mV", nominal = 1e-3);

    type ElectricCurrent = Modelica.Units.SI.ElectricCurrent(displayUnit = "meq/min", nominal = 9.64853399 * 10 ^ 4 / 1000 / 60);

    type Fraction = Real(final quantity = "Fraction", final unit = "1", displayUnit = "%", nominal = 1e-2);

    type HydraulicConductance = Real(final quantity = "HydraulicConductance", final unit = "m3/(Pa.s)", displayUnit = "l/(mmHg.min)", nominal = 1e-3 / (133.322387415 * 60), min = 0);

    type HydraulicResistance = Real(final quantity = "HydraulicConductance", final unit = "(Pa.s)/m3", displayUnit = "(mmHg.min)/l", nominal = 1e+3 * 133.322387415 * 60, min = 0);

    type HydraulicCompliance = Real(final quantity = "HydraulicCompliance", final unit = "m3/Pa", displayUnit = "ml/mmHg", nominal = 1e-6 / 133.322387415);

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
        Medium =   Air)                                                                      "External environment" annotation (
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
        Medium =   Air)                                                                      "External environment" annotation (
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
model modelECMORespiratoryVR_BloodGasesTransport_LungVentilatorSCMV2
 extends modelECMORespiratoryVR.BloodGasesTransport.LungVentilatorSCMV2;
  annotation(experiment(StopTime=10, __Dymola_Algorithm="Dassl"),uses(modelECMORespiratoryVR(version="1")));
end modelECMORespiratoryVR_BloodGasesTransport_LungVentilatorSCMV2;
