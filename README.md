# SBML files supported
Support is limited to SBML files meeting the following requirements:

1. SBML files of level 2, version 2 or 3 are supported. 
2. The file contains at least one reaction.
3. Each species must specify at least one of the values `InitialAmount` or `InitialConcentration`. If both are provided, `InitialAmount` takes precedence.
4. One or more compartments are supported, provided they do not appear in the kinetic law.
5. The presence of the <listOfParameters> tag is supported.
6. Only mass action kinetic laws are supported, and they must be represented in the file as follows:
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
#
| **`type` attribute** | **Example content** | **Description**                            |
| -------------------- | ------------------- | ------------------------------------------ |
| `integer`            | `5`                 | Integer number                             |
| `real`               | `3.14`              | Real number in decimal notation            |
| `real`               | `1.5e3`             | Scientific notation (exponential format)   |
| `rational`           | `3 4`               | Rational number as a fraction (3/4)        |
| `real`               | `2.5e+6`            | Scientific notation with positive exponent |
| `real`               | `1.0e-3`            | Scientific notation with negative exponent |
    

