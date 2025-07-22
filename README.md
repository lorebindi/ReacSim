# ReacSim
ReacSim is a chemical reaction simulator based on Gillespie Stochastic Simulation Algorithm (SSA), 
extended with support for events with and without delay. Additionally, it integrates some feature 
of LibRoadRunner to provide deterministic simulations using Ordinary Differential Equations (ODE).
## Run ReacSim Without Installation

To use **ReacSim** without installing it as a package, clone this repository and run the `main.py` script with Python using the following command-line arguments:

###  Command-Line Arguments

- `--filesbml <path_1> <path_2> ... <path_n>`  
  Specifies one or more SBML file paths to process. Both relative and absolute paths are supported.

- `--max-time-gillespie <time>`  
  Defines a timeout (in seconds) for the Gillespie simulation (optional)
  - If `time > 0`:
    - The Gillespie simulation runs with a timeout.
    - If the simulation completes within the timeout, the result is used.
    - If it exceeds the timeout, the process is terminated and the simulation continues using ODEs.
  - If `time == 0` or this option is omitted, both Gillespie and ODE simulations are executed without any timeout.
- `--t_max <max_value>` 
  Define max value for the time of simulation.


## SBML files supported
Support is limited to SBML files meeting the following requirements:

- SBML files of level 2, version 3 or 4 are supported. 
- The file must contains at least one reaction. Reversible reactions are supported.
- Each species must specify only the value `InitialAmount` otherwise it abort.
- One or more compartments are supported, provided they do not appear in the kinetic law.
- The presence of theTest/lotka_volterra_p_1.xml `<listOfParameters>` tag is supported only in the radix of the file.
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

## Gillespie SSA Extended with SBML Events

This custom implementation extends the standard Gillespie Stochastic Simulation Algorithm (SSA) to support **SBML events**.

In standard Gillespie simulation, system dynamics are entirely driven by reactions, which describe continuous and probabilistic transitions.  
In contrast, **SBML events** model discrete and instantaneous changes in the system (e.g., triggering a signal or instant concentration jump).  
Due to their fundamentally different nature, events and reactions must be handled in distinct phases of the simulation.


### Simulation Workflow Overview

At each iteration, the simulation checks whether the simulation time `t` exceeds the maximum time `t_max`.  
After computing the propensities, total rate `a_0`, and the time for the next reaction `tau`, the system proceeds as follows:

1. **Evaluation of previously scheduled delayed events**
2. **Detection of newly triggered events**
3. **Execution of valid events**
4. **Standard Gillespie reaction step** (only if no event is executed)


### 1. Evaluation of Scheduled Delayed Events

At every iteration, the simulator checks the list of delayed events scheduled in previous steps (`pending_event_delay`).  
For each event:

- It verifies whether the trigger condition still returns `true` at the scheduled `delay_time`:
  - if the trigger is now `false`, the event is considered invalid and is marked for removal;
  - if the trigger remains `true`, the simulator checks whether the event is scheduled in the interval `(t, t + τ]`;
    - if `delay_time` is the smallest found so far, it is stored in `min_delay_time` and the event is added to `to_eval_events`;
    - if equal to `min_delay_time`, the event is also added to `to_eval_events`.

Once this evaluation is complete:

- Events with invalid triggers are removed from `pending_event_delay`.
- If any delayed event must be executed now, simulation time `t` is advanced to `min_delay_time`.
- The list `to_eval_events` will contain all valid delayed events to be executed at this step.

If the event is invalid and the `useValuesFromTriggerTime` attribute was `true`, the corresponding cache (`event[VALUES_FROM_TRIGGER_TIME]`) is also cleared to avoid using outdated values.

### 2. Detection of Newly Triggered Events

The parser checks all model events and evaluates their `<trigger>` condition.  
An event is triggered **only if**:

- It was *false* in the previous step (`PREVIOUS = False`)
- It is *true* in the current step

This mechanism is controlled using the boolean flag `PREVIOUS`, which prevents repeated triggers of the same event.

Depending on the type of event:

- **Immediate events** (with no delay or zero delay) are directly added to `to_eval_events` to be executed during this step.
- **Delayed events** are stored in `pending_event_delay` and scheduled at time `t + delay_time`.  
  If `useValuesFromTriggerTime = true`, the values of involved species and parameters are stored in `event[VALUES_FROM_TRIGGER_TIME]`, to be restored when the event is executed.

### 3. Execution of Events

After event evaluation, if `to_eval_events` contains one or more events, they are executed before proceeding with the reaction logic.

For each event in the queue:

- **Trigger check**: The trigger condition is re-evaluated to ensure it still holds at the current time `t`.
- **State restoration (if required)**: If `useValuesFromTriggerTime = true`, the cached values (species and parameters) stored at trigger time are restored into `self.parser.species` and `self.parser.parameters`.
- **Event assignment execution**: Each `<eventAssignment>` is dynamically evaluated using `evaluate_expr`, and the resulting value is assigned to the corresponding variable.
  - If a computed value is negative, an exception is raised, and the event is not applied.

### 4. Reaction Selection (Standard SSA Step)

If no event is executed in the current step, the simulation proceeds with the standard Gillespie logic:

- A reaction is selected probabilistically based on current propensities.
- The system is updated accordingly.
- Time is incremented by `τ`.

