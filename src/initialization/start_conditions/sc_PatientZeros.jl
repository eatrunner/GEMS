export PatientZeros
export pathogen, ags
export initialize!


"""
    PatientZeros <: StartCondition

A `StartCondition` that infects a single individual in each of the given AGS (community identification numbers) at the beginning of the simulation.

# Fields
- `pathogen::String`: The pathogen with which the individual has to be infected
- `ags::Vector{Int64}`: A vector of AGS (community identification number) where the initial seeds should be planted

Example
```julia
# single individual infected at random in regions 13003000 and 13076033
# needs to be used with a population that contains these regions, e.g., "MV"
condition = PatientZeros(ags=[13003000, 13076033]) 
sim = Simulation(population = "MV", start_condition = condition)
```
"""
struct PatientZeros <: StartCondition
    pathogen::String
    teryt::Vector{Int64}

    function PatientZeros(;pathogen::String = "", teryt::Vector{Int64} = Int64[])
        isempty(teryt) && throw(ArgumentError("At least one ags must be provided!"))
        length(pathogen) > 0 && @warn "GEMS currently only supports single-pathogen simulations. Specifying a pathogen in PatientZeros will have no effect."
        # teryt.(teryt) # try casting to AGS (will throw error if invalid)
            
        return new(pathogen, teryt)
    end
end
Base.show(io::IO, cnd::PatientZeros) = write(io, "PatientZeros(AGS=$(join(cnd.ags, ",")))")


"""
    pathogen(patientzeros::PatientZeros)

Returns pathogen used to infect individuals at the beginning in this start condition.
"""
function pathogen(patientzeros::PatientZeros)
    return patientzeros.pathogen
end

"""
    ags(patientzeros::PatientZeros)::Vector{Int64}

Returns the vector of ags where intial seeds should be planted.
"""
function teryt(patientzeros::PatientZeros)::Vector{Int64}
    return patientzeros.teryt
end


"""
    initialize!(simulation::Simulation, condition::PatientZeros)

Initializes the simulation model, infecting a single individual in each of the regions provided by their AGS (community identification number).
"""
function initialize!(simulation::Simulation, condition::PatientZeros; seed_sample::Union{Int64,Nothing}=nothing)
    # create a new Xoshiro RNG for sampling, seeded from rng(simulation) if seed_sample is nothing, or from seed_sample otherwise
    rng_sample = isnothing(seed_sample) ? rng(simulation) : Xoshiro(seed_sample)
    # TODO handle pathogen selection
    # TODO handle multiple pathogens

    # number of individuals to infect
    to_infect = []
    for a in teryt(condition)
        # Get all individuals in households with the given ags
        inds = []
        for h in settings(simulation, Household)
            if a == h |> teryt
                push!(inds, individuals(h)...)
            end
        end
        if length(inds) == 0
            error("teryt($a) not in simulation or no individuals found in this region.")
        end
        # Sample one individual from the list of individuals
        append!(to_infect, gems_sample(rng_sample, inds, 200, replace=false))
    end
   
    # infect individuals
    for i in to_infect
        infect!(i, tick(simulation), pathogen(simulation), sim = simulation, rng = rng_sample)
        for (type, id) in settings(i, simulation)
            activate!(settings(simulation, type)[id])
        end
    end
end