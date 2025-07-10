# ReacSim
ReacSim is a chemical reaction simulator based on Gillespie Stochastic Simulation Algorithm (SSA), 
extended with support for events with and without delay. Additionally, it integrates some feature 
of LibRoadRunner to provide deterministic simulations using Ordinary Differential Equations (ODE).

## SBML files supported
Support is limited to SBML files meeting the following requirements:

- SBML files of level 2, version 3 or 4 are supported. 
- The file must contains at least one reaction. Reversible reactions are supported.
- Each species must specify only the value `InitialAmount` otherwise it abort.
- One or more compartments are supported, provided they do not appear in the kinetic law.
- The presence of the `<listOfParameters>` tag is supported only in the radix of the file.
- Only mass action kinetic laws are supported, and they must be represented in the file as follows:
```xml
<reaction id="R" reversible="false">
        <listOfReactants> ... </listOfReactants>
        <listOfProducts>  ... </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
            <times/>
            <ci> ... </ci>                  <!-- kinetic costant -->
            <apply>                         <!-- [<ci>]^<cn> -->
              <power/>
              <ci> ... </ci>
              <cn type="type"> ... </cn>
            </apply>
            <apply>                         <!-- [<ci>]^<cn> -->
              <power/>
              <ci> ... </ci>
              <cn type="type"> ... </cn>
            </apply>
            . 
            .
            .
          </apply>
          </math>
        </kineticLaw>
</reaction>
```

- Events with and without delay are supported.
  - The attribute `useValuesFromTriggerTime` is supported for L2V4 SBML files. In SBML L2V3 this attribute is not available but the documentation assumes the behaviour as if `useValuesFromTriggerTime = True`.    
  - Only addition and subtraction of constants and/or identifiers are supported in the `eventAssignment`'s expressions.

## Stochastic rate constant inference
ReacSim provide a way to infer the stochastic rate constant of a reaction starting from sperimental/simulated data. To enable this, it's important to follow these step:

- Provide an SBML file where the stochastic rate constant of the reaction appears in its kinekitc law of the reaction but is not listed in the `listOfParameters`.
- Place a csv file in the `Inference_of_stochastic_rate_constant` directory, named with the id of the constant in the SBML file. This file must contain:
  - one column for each molecular species involved in the reaction,
  - one column for time,
  - Several rows representing the discrete amounts of each molecular species at specific time points.

These steps automatically enable the stochastic rate constant inference during the validation of the kinetic law of that reaction.

The simulator also provides the `export_mean_species_counts_csv`, which generates a csv file containing simulated data for a specific reaction by performing several executions of the Gillespie SSA algorithm. 

### Inference process
The inference process works as follows:

- **Rate estimation**: At each time step, the simulator estimates the reaction rate based on species concentration changes, averaging across all species involved to reduce measurement noise.
- **Mass-action term calculation**: At each time step, the simulator computes the product of the reactants’ discrete amounts raised to their stoichiometric powers.
- **Linear regression**: Using a least-squares fit of the form `v = k * [A]^a * [B]^b ...`, the stochastic rate constant k is estimated.

## Gillespie SSA extended with SBML events 
