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
barn2cm2(b::Float64) = b * 1e-24


# Sim 
i = 1
nsim = 200

rMin = 2e-10
rMax = 1e-9

endNeutronCounts = zeros(nsim)
endActiveCounts  = zeros(nsim)
endOutsideCounts = zeros(nsim)

for sphereRadius in range(rMin, rMax, length = nsim)
    global i # counter for radius vectors
    println("[!] Radius $(sphereRadius)")

    # Experiment Parameters
    cycles = 15
    startingNeutrons = 500

    ## We assume an averaged energy for our neutrons, and the mean paths below 
    ##  correspond to that average energy. We take this average neutron energy 
    ##  to be thermal (in line with the goal of having a neutron regulator
    ##  inside fission reactors)
    ## Units are in Barns (= 10^(-24) cm^2)
    ## σ_sc : scattering cross section
    ## σ_cp : capture cross section
    ## σ_f  : fission cross section

    # Uranium-235 cross section averages for thermal neutrons (2200 m/s)
    σ_sc =  10.0
    σ_cp = 99.0
    σ_f = 583.0

    # Plutonium-239 cross section averages for thermal neutrons
    # σ_sc = 8.0
    # σ_cp = 269.0
    # σ_f = 748.0


    # Simulation

    ## Initialisation
    N = startingNeutrons # Neutron count
    V = (4/3) * π * sphereRadius^3 # We're in a sphere

    n = N / V           # neutron density

    σ_t = σ_sc + σ_cp + σ_f

    λ_sc = 1 / (σ_sc * n)
    λ_cp = 1 / (σ_cp * n)
    λ_f  = 1 / (σ_f * n)
    λ_t  = 1 / (barn2cm2(σ_t) * n) # total mean free path

    neutrons = [Neutron() for _ in 1:startingNeutrons]

    ## Data Collection
    insideNeutronCounts     = zeros(cycles)
    freeNeutronCounts       = zeros(cycles)
    outsideNeutronCounts    = zeros(cycles)
    absorbedNeutronCounts   = zeros(cycles)
    totalNeutronCounts      = zeros(cycles)

    ## Simulation

    for run in 1:cycles

        if run > 1
            outsideNeutronCounts[run]  = outsideNeutronCounts[run - 1]  # cumulative count
            absorbedNeutronCounts[run] = absorbedNeutronCounts[run - 1] # cumulative count
        end

        for j in 1:N
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
                outsideNeutronCounts[run] += 1
                continue
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
                absorbedNeutronCounts[run] += 1
            else 
                # Fission collision
                ## We just assume (n,2n) fission here

                # New neutron is released from fission, create it at this position
                newNeutron = Neutron(neutron.x, neutron.y, neutron.z)
                push!(neutrons, newNeutron)
                N += 1
            end
        end

        # Data Collection
        insideNeutronCounts[run] = N - outsideNeutronCounts[run]
        freeNeutronCounts[run] = N - outsideNeutronCounts[run] - absorbedNeutronCounts[run]
        totalNeutronCounts[run] = N
    end

    endNeutronCounts[i] = totalNeutronCounts[end]
    endActiveCounts[i]  = freeNeutronCounts[end]
    endOutsideCounts[i] = outsideNeutronCounts[end]
    i += 1
end




# Results Production
println("--== SIMULATION FINISHED ==--")
using Plots

x_range = range(rMin, rMax, length = nsim)
p = plot(x_range, endNeutronCounts, lc=:black, label = "Total")
plot!(x_range, endActiveCounts, lc=:green, label = "Active / Free")
plot!(x_range, endOutsideCounts, lc=:red, label = "Escaped")

xlabel!("Sphere Radius (cm)")
ylabel!("Neutron Count")
title!("U-235 @ N0=500, nSim=$(nsim)")

savefig(p, "results/radius-variation-uranium-semi-zoomed.pdf")