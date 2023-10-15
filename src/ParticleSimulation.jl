"""
This runs a particle simulation for `n` neutrons to generate histories
    and derive data from such


We assume isotropic scattering, an averaged neutron energy 
"""

# Data Structures and Helper Functions
mutable struct Neutron
    x::Float64
    y::Float64
    z::Float64

    escaped::Bool   # is this neutron outside the boundary?
    absorbed::Bool  # has this neutron been consumed by a nucleus? 

    # Default constructor
    Neutron() = new(0.0, 0.0, 0.0, false, false)
    Neutron(x::Float64, y::Float64, z::Float64) = new(x, y, z, false, false)
end

vec_length(x,y,z) = sqrt(x^2 + y^2 + z^2)

# Experiment Parameters
cycles = 30
startingNeutrons = 1000
maxNeutronsSupportable = 1e10

sphereRadius = 17.32 # cm radius

## We assume an averaged energy for our neutrons, and the mean paths below 
##  correspond to that average energy. We take this average neutron energy 
##  to be thermal (in line with the goal of having a neutron regulator
##  inside fission reactors)
## Units are in Barns (= 10^(-24) cm^2)
## σ_sc : scattering cross section
## σ_cp : capture cross section
## σ_f  : fission cross section

# Uranium-235 cross section averages for thermal neutrons (2200 m/s)
σ_sc =  10
σ_cp = 99
σ_f = 583

# Plutonium-239 cross section averages for thermal neutrons
# σ_sc = 8
# σ_cp = 269
# σ_f = 748


# Simulation

## Initialisation
N = startingNeutrons # Neutron count
V = sphereRadius

n = N / V           # neutron density

σ_t = σ_sc + σ_cp + σ_f

λ_sc = 1 / (σ_sc * n)
λ_cp = 1 / (σ_cp * n)
λ_f  = 1 / (σ_f * n)
λ_t  = 1 / (σ_t * n) # total mean free path

neutrons = [Neutron() for _ in 1:startingNeutrons]

## Data Collection


## Simulation
for run in 1:cycles
    for j in 1:N
        global N
        neutron = neutrons[j]

        if neutron.escaped || neutron.absorbed
            # Neutron is not in simulation anymore
            continue
        end

        # Determine travel direction
        L = -λ_t * log(rand())
        θ = asin(-1 + 2 * rand())
        ϕ = 2 * π * rand()

        # Calculate displacement from current (x,y,z)
        dx = L * cos(θ) * cos(ϕ)
        dy = L * cos(θ) * sin(ϕ)
        dz = L * sin(θ)

        # Update neutron position
        neutron.x = neutron.x + dx
        neutron.y = neutron.y + dy
        neutron.z = neutron.z + dz

        # Has the neutron escaped our medium?
        if vec_length(neutron.x, neutron.y, neutron.z) > sphereRadius
            neutron.escaped = true
        end

        # Random choice as to what collision event occurs
        r = rand()

        if r < (σ_sc / σ_t) 
            # Scattering collision
            ## Direction is chosen randomly at start of processing in each cycle,
            ##  and this is equivalent to a random direction, so we can continue here
            continue            
        elseif r < (σ_cp / σ_t)
            # Capture collision
            ## Neutron is absorbed safely by the nucleus
            ## NOTE: maybe some alpha / beta radiation or similar is emitted, 
            ##  or maybe that there is a delayed fission reaction, but here 
            ##  we ignore that and just assume our neutron disappears
            neutron.absorbed = true
        else 
            # Fission collision
            ## We just assume (n,2n) fission here

            # New neutron is released from fission, create it at this position
            newNeutron = Neutron(neutron.x, neutron.y, neutron.z)
            push!(neutrons, newNeutron)
            N += 1
        end

    end


end


# Results Production
println("--== SIMULATION FINISHED ==--")
println("Neutron Count: $(N)")
println("Outside      : $(length(filter(n -> n.escaped, neutrons)))")
println("Absorbed     : $(length(filter(n -> n.absorbed, neutrons)))")
println("Inside       : $(length(filter(n -> !(n.absorbed | n.escaped), neutrons)))")