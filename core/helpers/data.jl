
function get_timeline(xl::StateList, l::AbstractString)
    if l[1:8] == "microbes"
        test = false
        mark = 10
        while test == false
            if l[mark+1] == ']' 
                test = true
            else
                mark += 1
            end
        end
        s = Symbol(l[mark+3:end])
        tl = Float64[]
        for i = 1:length(xl)
            push!(tl, getfield(xl[i].microbes[eval(parse(l[10:mark]))], s))
        end
        return tl::Array{Float64, 1}
    elseif l[1:8] == "ecagents"
        test = false
        mark = 10
        while test == false
            if l[mark+1] == ']' 
                test = true
            else
                mark += 1
            end
        end
        s = Symbol(l[mark+3:end])
        tl = Float64[]
        for i = 1:length(xl)
            push!(tl, getfield(xl[i].ecagents[eval(parse(l[10:mark]))], s))
        end
        return tl::Array{Float64, 1}
    else
        s = Symbol(l)
        if typeof(getfield(xl[1],s)) == Float64
            tl = Float64[]
            for i = 1:length(xl)
                push!(tl, getfield(xl[i],s))
            end
            return tl::Array{Float64, 1}
        elseif typeof(getfield(xl[1],s)) == Array{Float64, 1}
            tl = getfield(xl[1], s)
            for i = 2:length(xl)
                tl = hcat(tl, getfield(xl[i], s))
            end
            return tl'::Array{Float64, 2}
        end
    end
end

function get_characteristic_time(xl::StateList, l::AbstractString, pcl::Array{Float64, 1}, g::Grid)
    vtl = get_timeline(xl, l)
    test = false
    for i = 1:length(vtl)-1
        if vtl[i] < vtl[i+1]
            test = true
        end
    end
    if test == false
        warn("characterictic times are not relevant: timeserie $l is not decreasing")
    end
    tl = Float64[]
    for pc in pcl
        i = 1
        while vtl[1] * (1- pc) < vtl[i] && i < length(vtl)
            i += 1
        end
        push!(tl, g.t[i])
    end
    return tl::Array{Float64, 1}
end

function save(xl::StateList, u::Context, p::Parameters; path="/", folder="")
    g = p.grid
    t = p.time
    mtl = p.microbes_traits
    ecl = p.ecagent_traits
    if path[end] != "/"
        path = path * "/"
    end

    
    if folder == ""
        folder = "results-"* Dates.format(now(), "yyyy-mm-dd\\THH:MM:SS")
        mkpath(path* folder)
    else
        mkpath(path* folder)
    end
    dir = path* folder

    # parameters
    writecsv("$dir/parameter_time.csv", ["unit" "min" "max"])
    file = open("$dir/parameter_time.csv", "a")
    writecsv(file, [t.label t.min t.max])
    close(file)
    
    writecsv("$dir/parameter_grid.csv", ["variable" "grid"])
    file = open("$dir/parameter_grid.csv", "a")
    writecsv(file, ["t" g.t'])
    writecsv(file, ["q" g.q'])
    close(file)
    
    writecsv("$dir/parameter_biochemical_classes.csv", ["class" "min" "max"])
    file = open("$dir/parameter_biochemical_classes.csv", "a")
    for bc in p.biochemical_classes
        writecsv(file, [bc.label bc.min bc.max])
    end
    close(file)
    
    writecsv("$dir/parameter_microbes.csv", ["microbe" "signature" "assimilation" "efficiency" "mortality"])
    file = open("$dir/parameter_microbes.csv", "a")
    for i=1:length(mtl)
        writecsv(file, [mtl[i].name mtl[i].signature_function mtl[i].assimilation mtl[i].efficiency mtl[i].mortality])
    end
    close(file)
    
    writecsv("$dir/parameter_ecagents.csv", ["ecagent" "kernel" "kernel_type" "activity"])
    file = open("$dir/parameter_ecagents.csv", "a")
    for i=1:length(ecl)
        writecsv(file, [ecl[i].name ecl[i].f_kernel ecl[i].kernel_type ecl[i].activity])
    end
    close(file)
    
    # context
    writecsv("$dir/context_input_dist.csv", ["t" "a"])
    file = open("$dir/context_input_dist.csv", "a")
    writecsv(file, u.input_dist)
    close(file)
    
    # state
    field_list_mic = ["uptake_flux", "assimilation_flux", "respiration_flux", "mortality_flux", "respiration", "living_carbon_mass", "dead_carbon_mass"]
    field_list_eca = ["action_rate", "transformation_flux"]
    variable_list = ["substrate_dist", "substrate", "substrate_classes", "substrate_classes_proportion", "respiration_flux", "respiration"]
    for i=1:length(xl[1].microbes)
        for field in field_list_mic
            push!(variable_list, "microbes[$i].$field")
        end
    end
    for i=1:length(xl[1].ecagents)
        for field in field_list_eca
            push!(variable_list, "ecagents[$i].$field")
        end
    end
    
    for variable in variable_list
        data = get_timeline(xl, variable)
        # create header
        writecsv("$dir/state_$variable.csv", ["time" variable])
        # add data
        file = open("$dir/state_$variable.csv", "a")
        writecsv(file, [g.t data])
        close(file)
    end
end

function load(pth::AbstractString)
    data = readcsv(pth* "/parameter_biochemical_classes.csv")
    bcc = BiochemicalClass[]
    for i = 2:size(data)[1]
        push!(bcc, BiochemicalClass(data[i,1], convert(Float64, data[i,2]), convert(Float64, data[i,3])))
    end
    
    data = readcsv(pth* "/parameter_time.csv")
    time = Time(data[2,1], convert(Float64, data[2,2]), convert(Float64, data[2,3]))
    
    data = readcsv(pth* "/parameter_grid.csv")
    g = Grid(bcc, data[3,3] - data[3,2], time, data[2,3] - data[2,2])

    data = readcsv(pth* "/parameter_ecagents.csv")
    etl = ECAgentTraits[]
    for i=2:size(data)[1]
        et = ECAgentTraits(data[i,1], data[i,2], data[i,3], data[i,4], g)
        push!(etl, et)
    end

    data = readcsv(pth* "/parameter_microbes.csv")
    mtl = MicrobeTraits[]
    for i=2:size(data)[1]
        mt = MicrobeTraits(data[i,1], data[i,2], data[i,3], data[i,4], data[i,5], g)
        push!(mtl, mt)
    end

    p = Parameters(bcc, time, g, mtl, etl) 
    
    data = readcsv(pth* "/context_input_dist.csv")[2:end, 2:end]
    for i=1:size(data)[1]
        for j=1:size(data)[2]
            data[i,j] = convert(Float64, data[i,j])
        end
    end
    u = Context(data)
    xl = State[]

    substrate_temp = convert(Array{Float64, 1}, readcsv(pth* "/state_substrate.csv")[2:end, 2])
    substrate_dist_temp = convert(Array{Float64, 2}, readcsv(pth* "/state_substrate_dist.csv")[2:end, 2:end])
    substrate_classes_temp = convert(Array{Float64, 2}, readcsv(pth* "/state_substrate_classes.csv")[2:end, 2:end])
    substrate_classes_proportion_temp = convert(Array{Float64, 2}, readcsv(pth* "/state_substrate_classes_proportion.csv")[2:end, 2:end])
    respiration_flux_temp = convert(Array{Float64, 1}, readcsv(pth* "/state_respiration_flux.csv")[2:end, 2])
    respiration_temp = convert(Array{Float64, 1}, readcsv(pth* "/state_respiration.csv")[2:end, 2])
    
    nmic = size(readcsv(pth* "/parameter_microbes.csv"))[1]- 1
    neca = size(readcsv(pth* "/parameter_ecagents.csv"))[1]- 1
    mic_uptake_flux_temp = Array{Float64, 1}[]
    mic_assimilation_flux_temp = Array{Float64, 1}[]
    mic_respiration_flux_temp = Array{Float64, 1}[]
    mic_mortality_flux_temp = Array{Float64, 1}[]
    mic_respiration_temp = Array{Float64, 1}[]
    mic_living_carbon_mass_temp = Array{Float64, 1}[]
    mic_dead_carbon_mass_temp = Array{Float64, 1}[]
    for j = 1:nmic
        push!(mic_uptake_flux_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].uptake_flux.csv")[2:end, 2]))
        push!(mic_assimilation_flux_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].assimilation_flux.csv")[2:end, 2]))
        push!(mic_respiration_flux_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].respiration_flux.csv")[2:end, 2]))
        push!(mic_mortality_flux_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].mortality_flux.csv")[2:end, 2]))
        push!(mic_respiration_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].respiration.csv")[2:end, 2]))
        push!(mic_living_carbon_mass_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].living_carbon_mass.csv")[2:end, 2]))
        push!(mic_dead_carbon_mass_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_microbes[$j].dead_carbon_mass.csv")[2:end, 2]))
    end
    eca_action_rate_temp = Array{Float64, 1}[]
    eca_transformation_flux_temp = Array{Float64, 1}[]
    for j = 1:neca
        push!(eca_action_rate_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_ecagents[$j].action_rate.csv")[2:end, 2]))
        push!(eca_transformation_flux_temp, convert(Array{Float64, 1}, readcsv(pth* "/state_ecagents[$j].transformation_flux.csv")[2:end, 2]))
    end
    
    for i = 1:length(g.t)
        mxl = MicrobeState[]
        for j = 1:nmic
            mic_uptake_flux = mic_uptake_flux_temp[j][i]
            mic_assimilation_flux = mic_assimilation_flux_temp[j][i]
            mic_respiration_flux = mic_respiration_flux_temp[j][i]
            mic_mortality_flux = mic_mortality_flux_temp[j][i]
            mic_respiration = mic_respiration_temp[j][i]
            mic_living_carbon_mass = mic_living_carbon_mass_temp[j][i]
            mic_dead_carbon_mass = mic_dead_carbon_mass_temp[j][i]
            mx = MicrobeState(mic_uptake_flux, mic_assimilation_flux, mic_respiration_flux, mic_mortality_flux, mic_respiration, mic_living_carbon_mass, mic_dead_carbon_mass)
            push!(mxl, mx)
        end
        exl = ECAgentState[]
        for j = 1:neca
            eca_action_rate = eca_action_rate_temp[j][i]
            eca_transformation_flux = eca_transformation_flux_temp[j][i]
            ex = ECAgentState(eca_action_rate, eca_transformation_flux)
            push!(exl, ex)
        end
        substrate = substrate_temp[i]
        substrate_dist = substrate_dist_temp[i, :]
        substrate_classes = substrate_classes_temp[i, :]
        substrate_classes_proportion = substrate_classes_proportion_temp[i, :]
        respiration_flux = respiration_flux_temp[i]
        respiration = respiration_temp[i]
        x = State(substrate_dist, substrate, substrate_classes, substrate_classes_proportion, respiration_flux, respiration, mxl, exl)
        push!(xl, x)
    end

    return p::Parameters, u::Context, xl::StateList
end
