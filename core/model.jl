
function get_index(array::Array, subarray::Array)
    index = Int[]
    for s in subarray
        i = 1
        while array[i] < s && i < length(array)
            i += 1 
        end
        push!(index, i)
    end
    return index::Array{Int, 1}
end

type Interval
    min::Float64
    max::Float64
    norm::Float64
    function Interval(min::Float64, max::Float64)
        if max-min<0
            error("invalid LabeledInterval bounds")
        end
        new(min, max, max- min)
    end
    function Interval(a::Array{Float64,1})
        if a[2]-a[1]<0
            error("invalid LabeledInterval bounds")
        end
        new(a[1], a[2], a[2]- a[1])
    end
end

const IntervalList = Array{Interval, 1}

const Domain = Interval

const DomainList = IntervalList

type LabeledInterval
    label::AbstractString
    min::Float64
    max::Float64
    norm::Float64
    function LabeledInterval(label::AbstractString, min::Float64, max::Float64)
        if max-min<0
            error("invalid LabeledInterval bounds")
        end
        new(label, min, max, max-min)
    end
end

const LabeledIntervalList = Array{LabeledInterval, 1}

function get_interval(li::LabeledInterval)
    return Interval(li.min, li.max)::Interval
end

function label_interval(l::AbstractString, i::Interval)
    return LabeledInterval(l, i.min, i.max)
end

function get_dictionnary(lil::LabeledIntervalList)
    d = Dict{AbstractString, Interval}()
    for li in lil
        d[li.label] = get_interval(li)
    end
    return d
end

const BiochemicalClass = LabeledInterval

const BiochemicalClasses = LabeledIntervalList

const Time = LabeledInterval

const Timeline = Array{Float64, 1}

type Grid
    qmin::Float64
    qmax::Float64
    dq::Float64
    q::Array{Float64, 1}
    integrate::Function
    integrate_class::Function
    integrate_interval::Function
    tmin::Float64
    tmax::Float64
    dt::Float64
    t::Array{Float64, 1}
    function Grid(bcc::BiochemicalClasses, 
                  dq::Float64,
                  time::Time,
                  dt::Float64)
        qmin = min([bcc[i].min for i=1:length(bcc)]...)
        qmax = max([bcc[i].max for i=1:length(bcc)]...)
        q = qmin:dq:qmax
        t = collect(time.min:dt:time.max)
        function integrate(array::Array{Float64, 1})
            return sum(array[1:end-1])* dq
        end
        function integrate_class(rho::Array{Float64, 1}, 
                                 bcc::BiochemicalClass)
            return integrate([rho[i] for i in get_index(collect(q), collect(bcc.min:dq:bcc.max))])
        end
        function integrate_interval(rho::Array{Float64, 1}, 
                                    i::Interval)
            return integrate([rho[j] for j in get_index(collect(q), collect(i.min:dq:i.max))])
        end
        new(qmin, qmax, dq, q, integrate, integrate_class, integrate_interval, time.min, time.max, dt, t)
    end
end

type ECAgentTraits
    name::AbstractString
    f_kernel::AbstractString
    kernel::Array{Float64, 2}
    kernel_type::AbstractString
    activity::AbstractString
    activity_function::Function
    function ECAgentTraits(name::AbstractString,
                           f_kernel::AbstractString,
                           kernel_type::AbstractString,
                           activity::AbstractString,
                           g::Grid)
        function_kernel = eval(parse(f_kernel))
        activity_function = eval(parse(activity))
        if kernel_type == "regular" || kernel_type == "translation"
            kernel = zeros(length(g.q),length(g.q))
            for k=1:length(g.q)
                kernel[k,:] = [function_kernel(g.q[k], g.q[j]) for j=1:length(g.q)]
            end
        elseif kernel_type == "integral"
            kernel = zeros(length(g.q),length(g.q))
            for k=1:length(g.q)
                for j = 1:length(g.q)
                    if j==1
                        kernel[k,j] = 0.5 * Base.invokelatest(function_kernel, g.q[k], g.q[j], g)
                    elseif j == length(g.q)
                        kernel[k,j] = 0.5 * Base.invokelatest(function_kernel, g.q[k], g.q[j-1], g)
                    else
                        kernel[k,j] = 0.5 * (Base.invokelatest(function_kernel, g.q[k], g.q[j], g) +
                                             Base.invokelatest(function_kernel, g.q[k], g.q[j-1], g))
                    end
                end
            end
        end
        new(name, f_kernel, kernel, kernel_type, activity, activity_function)
    end
end

const ECAgentTraitsList = Array{ECAgentTraits, 1}

type MicrobeTraits
    name::AbstractString
    signature_function::AbstractString
    signature::Array{Float64, 1}
    assimilation::AbstractString
    assimilation_function::Function
    efficiency::AbstractString
    efficiency_function::Function
    mortality::AbstractString
    mortality_function::Function
    function MicrobeTraits(name::AbstractString,
                           signature_function::AbstractString,
                           assimilation::AbstractString,
                           efficiency::AbstractString,
                           mortality::AbstractString,
                           g::Grid)
        new(name, signature_function, [Base.invokelatest(eval(parse(signature_function)), g.q[i]) for i=1:length(g.q)], assimilation, eval(parse(assimilation)), efficiency, eval(parse(efficiency)), mortality, eval(parse(mortality)))
    end
end

const MicrobeTraitsList = Array{MicrobeTraits, 1}

type Parameters
    biochemical_classes::BiochemicalClasses
    time::Time
    grid::Grid
    microbes_traits::MicrobeTraitsList
    ecagent_traits::ECAgentTraitsList
end

type MicrobeState
    uptake_flux::Float64
    assimilation_flux::Float64
    respiration_flux::Float64
    mortality_flux::Float64
    respiration::Float64
    living_carbon_mass::Float64
    dead_carbon_mass::Float64
end

const MicrobeStateList = Array{MicrobeState, 1}

type ECAgentState
    action_rate::Float64
    transformation_flux::Float64
end

const ECAgentStateList = Array{ECAgentState, 1}

type State
    substrate_dist::Array{Float64, 1}
    substrate::Float64
    substrate_classes::Array{Float64, 1}
    substrate_classes_proportion::Array{Float64, 1}
    respiration_flux::Float64
    respiration::Float64
    microbes::MicrobeStateList
    ecagents::ECAgentStateList
    function State(substrate_dist::Array{Float64, 1},
                   respiration_flux::Float64,
                   respiration::Float64,
                   microbes::MicrobeStateList, 
                   ecagents::ECAgentStateList,
                   p::Parameters)
        substrate = p.grid.integrate(substrate_dist)
        substrate_classes = zeros(length(p.biochemical_classes))
        substrate_classes_proportion = zeros(length(p.biochemical_classes))
        for i = 1:length(p.biochemical_classes)
            substrate_classes[i] = p.grid.integrate_class(substrate_dist, p.biochemical_classes[i])
            substrate_classes_proportion[i] = substrate_classes[i]/ substrate
        end
        new(substrate_dist, substrate, substrate_classes, substrate_classes_proportion, respiration_flux, respiration, microbes, ecagents)
    end
    function State(substrate_dist::Array{Float64, 1}, substrate::Float64, substrate_classes::Array{Float64, 1}, substrate_classes_proportion::Array{Float64, 1}, respiration_flux::Float64, respiration::Float64, microbes::MicrobeStateList, ecagents::ECAgentStateList)
       new(substrate_dist, substrate, substrate_classes, substrate_classes_proportion, respiration_flux, respiration, microbes, ecagents)
    end
    function State()
       new()
    end
end

const StateList = Array{State, 1}

type Context
    input_dist::Array{Float64, 2}
end

function nextstate_explicit(xn::StateList, 
                            tstep::Int, 
                            p::Parameters, 
                            u::Context)
    t = p.grid.t[tstep]
    la = length(p.grid.q)
    lt = length(p.grid.t)
    nm = length(p.microbes_traits)
    ne = length(p.ecagent_traits)

    sub = xn[end].substrate_dist + p.grid.dt* u.input_dist[:,tstep]

    mic_ns_list = MicrobeState[]
    respiration_flux_t = 0.
    dsub_mic = zeros(sub)
    for i = 1:nm
        efficiency = [p.microbes_traits[i].efficiency_function(xn, u, p.grid.q[j], t) for j=1:la]
        uptake = [p.microbes_traits[i].assimilation_function(xn, u, p.grid.q[j], t) for j=1:la].* sub
        uptake_flux = p.grid.integrate(uptake)
        assimilation_flux = p.grid.integrate(uptake.* efficiency)
        respiration_flux = p.grid.integrate(uptake.* (1- efficiency))
        mortality_flux = p.microbes_traits[i].mortality_function(xn, u, t)* xn[end].microbes[i].living_carbon_mass

        living_carbon_mass = xn[end].microbes[i].living_carbon_mass+ p.grid.dt* (assimilation_flux- mortality_flux)
        dead_carbon_mass = xn[end].microbes[i].dead_carbon_mass+ p.grid.dt* mortality_flux
        respiration = xn[end].microbes[i].respiration+ p.grid.dt* respiration_flux

        mic_ns = MicrobeState(uptake_flux, assimilation_flux, respiration_flux, mortality_flux, respiration, living_carbon_mass, dead_carbon_mass)
        push!(mic_ns_list, mic_ns)
        respiration_flux_t += respiration_flux
        dsub_mic += p.grid.dt* (mortality_flux* p.microbes_traits[i].signature- uptake)
    end
    dsub = dsub_mic
    respiration_t = xn[end].respiration + p.grid.dt * respiration_flux_t

    eca_ns_list = ECAgentState[]
    mat_enz = zeros(la,la)
    for i = 1:ne
        action_rate = [p.ecagent_traits[i].activity_function(xn, u, p.grid.q[j], t) for j=1:la]
        action = action_rate .* sub
        kernel = p.ecagent_traits[i].kernel
        mat_enz -= p.grid.dt * diagm(action_rate)
        if p.ecagent_traits[i].kernel_type == "regular"
            kernel = p.grid.dt * p.grid.dq * kernel
            for n = 1:la-1
                kernel[:,n] = kernel[:,n] * action_rate[n]
            end
        elseif p.ecagent_traits[i].kernel_type == "integral"
            kernel = p.grid.dt * kernel
            for n = 1:la
                kernel[:,n] = kernel[:,n] * action_rate[n]
            end
        elseif p.ecagent_traits[i].kernel_type == "translation"
            kernel = p.grid.dt * kernel
            for n = 1:la-1
                kernel[:,n] = kernel[:,n] * action_rate[n]
            end
        end
        mat_enz += kernel
        
        eca_ns = ECAgentState(p.grid.integrate(action_rate), p.grid.integrate(action))
        push!(eca_ns_list, eca_ns)
    end
    dsub_enz = mat_enz * sub
    dsub += dsub_enz

    sub_ns = sub + dsub 
    return State(sub_ns, respiration_flux_t, respiration_t, mic_ns_list, eca_ns_list, p)::State
end

function simulate_explicit(x0::State, 
                           p::Parameters, 
                           u::Context,
                           disp::Bool)
    xnl = State[];
    push!(xnl, x0)
    for inc = 2:length(p.grid.t)
        if disp == true
            println(p.grid.t[inc])
        end
        new = nextstate_explicit(xnl, inc, p, u)
        push!(xnl, new)
    end
    return xnl::StateList
end

include("helpers/functions.jl")
include("helpers/data.jl")
println("--------------- core loaded ---------------")
