include("../core/model.jl")

# framework
bcc = BiochemicalClasses([BiochemicalClass("Microbial sugar", 6., 9.),                      
                          BiochemicalClass("Plant sugar", 9., 12.)]) 

t = Time("days", 0., 300.) 
g = Grid(bcc, 0.01, t, 0.1) 

# traits
alpha = 5.
tau_0 = 1.8

et = ECAgentTraits("enzyme",                   
                   "(q, qp, g::Grid) -> integral_kernel_q_increase(Domain(10., 12.), q, qp, $alpha, g)",
                   "integral",                   
                   "(xl::StateList, u::Context, q, t) -> uptake_mic(xl, Domain(10., 12.), $tau_0, 1, q)",
                   g)
etl = ECAgentTraits[et]

u_0 = 5.
m_0 = .02
e_0 = .4

mt = MicrobeTraits("microbe",
                   "q -> gaussian(Domain(6., 9.), 7.5, .1, q)",
                   "(xl::StateList, u::Context, q, t) -> uptake_mic(xl, Domain(11.6, 12.), $u_0, 1, q)",    
                   "(xl::StateList, u::Context, q, t) -> $e_0",
                   "(xl::StateList, u::Context, t) -> $m_0",
                   g) 
mtl = MicrobeTraits[mt] 

p = Parameters(bcc, t, g, mtl, etl) 

# initial state 
ratio_mic_tot = 0.05

sub_0 = (1- ratio_mic_tot)* [gaussian(get_interval(bcc[2]), 10.5, .1, p.grid.q[i]) for i=1:length(p.grid.q)] 

m_uptake_flux_0 = 0.
m_assimilation_flux_0 = 0.
m_respiration_flux_0 = 0.
m_mortality_flux_0 = 0.  
m_respiration_0 = 0. 
m_living_carbon_mass_0 = ratio_mic_tot
m_dead_carbon_mass_0 = 0.
ms_0 = MicrobeState(m_uptake_flux_0, 
                    m_assimilation_flux_0,
                    m_respiration_flux_0, 
                    m_mortality_flux_0, 
                    m_respiration_0, 
                    m_living_carbon_mass_0, 
                    m_dead_carbon_mass_0) 
msl_0 = MicrobeState[ms_0] 

es_0 = ECAgentState(0., 0.)
esl_0 = ECAgentState[es_0]

x0 = State(sub_0, 0., 0., msl_0, esl_0, p) 

# context 
input = zeros(length(p.grid.q), length(p.grid.t)) 
u = Context(input) 

# simulation 
@time xl = simulate_explicit(x0, p, u, false) 
whos(r"xl") 
save(xl, u, p; path = "results", folder="example")
