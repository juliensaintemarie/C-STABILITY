# Model of organic matter transformation by microbial and extra-cellular activities according to accesibility

## Software requirements

 - Install Julia 1.0- or above  https://julialang.org/
 - Install library LinearAlgebra, DelimitedFiles, Dates and InteractiveUtils

## How to simulate?

### Main principle of a simulation

Model variables are stored for each simulation time in a Julia Struct called State (noted x).

The results of a simulation are store in a StateList (noted xl) collecting the states for all simulation times.

The model funtionning is given by **xl = simulate(x0, p, c)** where
	
 1. *x0* is the *initial state* of the system,
 2. *p* are the simulation *parameters*, and
 3. *c* is the simulation *context* (e.g. environmental variables over the time).

#### Step 1 - define the parameters p

The definition of microbes and extra-cellular agents traits requires the use of functions.
Several functions are defined in "core/helpers/functions.jl".

#### Step 2 - define the initial state x0


#### Step 3 - define the context c


 
### Structure and Types

The main structures used to perform a simulation are:
 - Parameters - p
 - Context - c
 - State - x

Elementary Julia Types (Int, Float64, AbstractString, Bool, Function, ...) are not detailed here.

| Category 		 | New Struct   	  	        | Abr.	| Elements     	   	    	      	      	| Type or Struct					|
|------------------------|------------------------------|-------|-----------------------------------------------|---------------------------------------|
| Basic types		 | Interval			| i	| min						| Float64				|
|			 |				| -	| max						| Float64				|
|			 |				| -	| norm						| Float64				|
|			 | IntervalList  		| il	| 						| Array{Interval, 1}			|
|			 | LabeledInterval		| li	| label						| AbstractString			|
|			 | 				| -    	| min						| Float64				|
|			 | 				| -	| max						| Float64				|
|			 |				| -	| norm						| Float64				|
|			 | LabeledIntervalList 		| lil	| 						| Array{LabeledInterval, 1}		|
|			 | Domain 			| d	| 						| Interval		 		|
|			 | DomainList 			| dl	| 						| IntervalList				|
|			 | Timeline 			| tl	|						| Array{Float64, 1}			|
| Parameters		 | BiochemicalClass		| bc	|						| LabeledInterval			|
|			 | BiochemicalClasses		| bcl	| 						| LabeledIntervalList			|
|			 | Time 			| t	| 						| LabeledInterval			|
| 			 | Grid				| g	| qmin						| Float64				|
|			 |				| -	| qmax						| Float64				|
|			 |				| -	| dq						| Float64				|
|			 |				| -	| q						| Array{Float64, 1}			|
|			 |				| -	| integrate					| Function				|
|			 |				| -	| integrate_class				| Function				|
|			 |				| -	| integrate_interval				| Function				|
|			 |				| -	| tmin						| Float64				|
|			 |				| -	| tmax						| Float64				|
|			 |				| -	| dt						| Float64				|
|			 |				| -	| t						| Array{Float64, 1}			|
|			 | MicrobeTraits		| mt	| name						| AbstractString			|
|			 | 				| -	| signature					| Array{Float64, 1}			|
|			 | 				| -	| assimilation					| AbstractString			|
|			 | 				| -	| assimilation_function				| Function				|
|			 | 				| -	| efficiency					| AbstractString			|
|			 | 				| -	| efficiency_function				| Function				|
|			 | 				| -	| mortality					| AbstractString			|
|			 | 				| -	| mortality_function				| Function				|
|			 | MicrobeTraitsList 		| mtl 	| 						| Array{MicrobeTraits, 1}		|
|			 | ECAgentTraits		| et 	| name						| AbstractString			|
|			 | 				| -	| f_kernel					| AbstractString			|
|			 | 				| -	| kernel					| Array{Float64, 2}			|
|			 | 				| -	| kernel_type					| AbstractString			|
|			 | 				| -	| activity					| AbstractString			|
|			 | 				| -	| activity_function				| Function				|
|			 | ECAgentTraitsList 		| etl 	| 						| Array{ECAgentTraits, 1}		|
|			 | Parameters			| p 	| biochemical_classes  				| BiochemicalClasses			|
|			 | 				| -	| time						| Time					|
|			 |				| -	| grid						| Grid					|
|			 | 				| -	| microbes_traits				| MicrobeTraitsList			|
|			 | 				| -	| ecagent_traits				| ECAgentTraitsList			|
| Context		 | Context			| c	| input_dist					| Array{Float64, 2}			|
| State 		 | MicrobeState			| mx	| uptake_flux					| Float64				|
|			 |    				| -	| assimilation_flux				| Float64				|
|			 |    				| -	| respiration_flux				| Float64				|
|			 |    				| -	| mortality_flux				| Float64				|
|			 |    				| -	| respiration					| Float64				|
|			 |    				| -	| living_carbon_mass				| Float64				|
|			 |    				| -	| dead_carbon_mass				| Float64				|
|			 | MicrobeStateList 		| mxl	| 						| Array{MicrobeState, 1}		|
|			 | ECAgentState			| ex 	| action_rate					| Float64	      			|
|			 |    				| -	| transformation_flux				| Float64				|
|			 | ECAgentStateList 		| exl 	| 						| Array{ECAgentState, 1}		|
|			 | State			| x	| substrate_dist				| Array{Float64, 1}			|
|			 |				| -	| substrate					| Float64		 		|
|			 |				| -	| substrate_classes				| Array{Float64, 1}			|
|			 |				| -	| substrate_classes_proportion			| Array{Float64, 1}			|
|			 |				| -	| respiration_flux				| Float64	    			|
|			 |				| -	| respiration					| Float64				|
|			 |				| -	| microbes					| MicrobeStateList			|
|			 |				| -	| ecagents					| ECAgentStateList			|
|			 | StateList			| xl	| 						| Array{State, 1}			|

