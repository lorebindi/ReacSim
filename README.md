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

Please note that the `type`:

| **`type` attribute** | **Example content** | **Description**                            |
| -------------------- | ------------------- | ------------------------------------------ |
| `integer`            | `5`                 | Integer number                             |
| `real`               | `3.14`              | Real number in decimal notation            |
| `real`               | `1.5e3`             | Scientific notation (exponential format)   |
| `rational`           | `3 4`               | Rational number as a fraction (3/4)        |
| `real`               | `2.5e+6`            | Scientific notation with positive exponent |
| `real`               | `1.0e-3`            | Scientific notation with negative exponent |

- Events with and without delay are supported.
  - The attribute `useValuesFromTriggerTime` is supported for L2V3 SBML files.
  - Only addition and subtraction of constants and/or identifiers are supported.

## Stochastic rate constant inference

## Gillespie SSA extended with SBML events 
