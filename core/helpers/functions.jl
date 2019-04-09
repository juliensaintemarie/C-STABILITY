#### general functions

function indicatrice(min::Float64, max::Float64, x::Float64)
    out = 0.
    if min <= x && x < max
        out = 1.
    end
    return out::Float64
end

function indicatrice(i::Interval, x::Float64)
    return indicatrice(i.min, i.max, x)::Float64
end

function linear(i::Interval, rate::Float64, x::Float64)
    return indicatrice(i, x)* rate::Float64
end

function linear(min::Float64, max::Float64, rate::Float64, x::Float64)
    return linear(Interval(min, max), rate, x)::Float64
end

function gaussian(i::Interval, x0::Float64, sigma::Float64, x::Float64) 
    return indicatrice(i, x)/ sqrt(2.*pi)/sigma * exp(-1./2.*((x-x0)/sigma)^2)::Float64
end

#### uptakes

function uptake_mic(xl::StateList, d::Domain, rate::Float64, mic::Int, q::Float64)
    return linear(d, rate, q)* xl[end].microbes[mic].living_carbon_mass::Float64
end

function uptake_flux(xl::StateList, d::Domain, rate::Float64, eca::Int, q::Float64)
    return linear(d, rate, q)* xl[end].ecagents[eca].transformation_flux::Float64
end

function uptake_michaelis_mic(xl::StateList, d::Domain, rate::Float64, K::Float64, mic::Int, q::Float64)
    mass = xl[end].microbes[mic].living_carbon_mass
    return linear(d, rate, q)* mass/ (K+ mass)::Float64
end

#### competition

function competition(xl::StateList, mic_ref::Int, mic::Int, lambda::Float64)
    lambda = 1.
    c = 1.- 1./ (1.+ 9. * exp(-lambda * xl[end].microbes[mic].living_carbon_mass/ xl[end].microbes[mic_ref].living_carbon_mass))
    return c::Float64
end

#### extra-cellular agent kernel

function kernel_uniform(d::Domain, q::Float64, dp::Domain, qp::Float64)
    return indicatrice(d, q)* indicatrice(dp, qp)/ f.norm::Float64
end

function integral_kernel_uniform(d::Domain, q::Float64, dp::Domain, qp::Float64, g::Grid)
    return indicatrice(d, q)* indicatrice(dp.min, dp.max - g.dq, qp)* g.dq/ d.norm::Float64
end

function translation_kernel(ht::Float64, qi::Float64, d::Domain, qj::Float64, g::Grid)
    out = 0.
    if isapprox(qi, qj+ht)
        out = 1.
    end
    return out::Float64
end

function integral_kernel_translation(ht::Float64, qi::Float64, d::Domain, qj::Float64, g::Grid)
    out = 0.
    if isapprox(qi, qj+ht)  && qi < d.max+ ht
        out = 1.
    end
    return out::Float64
end

function kernel_q_increase(d::Domain, q::Float64, qp::Float64, alpha::Float64)
    min = d.min
    max = d.max
    if min<qp && qp<=q && q<max
        if max == qp
            out = 0.
        else
            out = (alpha+1)* (max-q)^alpha/ (max-qp)^(alpha+1)
        end
    else
        out = 0.
    end
    return out::Foat64
end

function integral_kernel_q_increase(d::Domain, qj::Float64, qi::Float64, alpha::Float64, g::Grid)
    min = d.min
    max = d.max
    qiplus1 = qi + g.dq
    if min<qi && qiplus1<=qj && qj<max
        out = (alpha+1) / alpha * (max-qj)^alpha * (1/(max-qiplus1)^alpha - 1/(max-qi)^(alpha))
    else
        out = 0.
    end
    return out::Float64
end

function kernel_q_increase_stabilized(min::Float64, max::Float64, q::Float64, qp::Float64, alpha::Float64, c::Float64, dq::Float64)
    if qp >= min && qp < c && q >= qp && q <= max
        out = (alpha+1)* (max-q)^alpha/ (max-qp)^(alpha+1)
    elseif qp >= c && qp <= max && q == qp
        out = 1/dq
    else
        out = 0.
    end
    return out::Float64
end
