## Required packages for this script
using Random, Distributions, Plots

# Set random seed for reproducibility
Random.seed!(123)

# Number of intermediate states
number_of_intermediate_states = 1;

# Define rate parameters
stem_cell_death_rate = 0.001
stem_cell_differentiation_rate = 0.2

differentiated_cell_death_rate = 0.01
differentiated_cell_dedifferentiation_rate = 0.05

# Function for the expected extinction time for a chains starting in state S
function no_intermediate_states_expected_extinction_time_from_S(stem_cell_death_rate,stem_cell_differentiation_rate,differentiated_cell_death_rate,differentiated_cell_dedifferentiation_rate)
    (differentiated_cell_dedifferentiation_rate + differentiated_cell_death_rate + stem_cell_differentiation_rate)/(stem_cell_differentiation_rate*differentiated_cell_death_rate + stem_cell_death_rate*(differentiated_cell_dedifferentiation_rate + differentiated_cell_death_rate))
end

## One intermediate state
# Additional parameters (compared to the no intermediate state case)
intermediate_differentiation_rate = 0.25
intermediate_dedifferentiation_rate = 0.1
intermediate_death_rate = 0.01

# Function for the expected extinction time for a chains starting in state S
function one_intermediate_states_expected_extinction_time_from_S(stem_cell_death_rate, stem_cell_differentiation_rate,
    intermediate_death_rate, intermediate_differentiation_rate, intermediate_dedifferentiation_rate,
    differentiated_cell_death_rate, differentiated_cell_dedifferentiation_rate)
    det_F = (differentiated_cell_dedifferentiation_rate + differentiated_cell_death_rate)*
    (stem_cell_death_rate*(intermediate_dedifferentiation_rate + intermediate_death_rate) +
    stem_cell_differentiation_rate*intermediate_death_rate) +
    intermediate_differentiation_rate*differentiated_cell_death_rate*
    (stem_cell_differentiation_rate + stem_cell_death_rate)
    return (1/det_F)*((differentiated_cell_dedifferentiation_rate + differentiated_cell_death_rate)*
    (intermediate_dedifferentiation_rate + intermediate_death_rate + stem_cell_differentiation_rate) + 
    intermediate_differentiation_rate*(stem_cell_differentiation_rate + differentiated_cell_death_rate))
end

# Calculate expected time to extinction for current parameters
# expected_extinction_time = no_intermediate_states_expected_extinction_time_from_S(stem_cell_death_rate,stem_cell_differentiation_rate,differentiated_cell_death_rate,differentiated_cell_dedifferentiation_rate)
expected_extinction_time = one_intermediate_states_expected_extinction_time_from_S(stem_cell_death_rate, stem_cell_differentiation_rate,
intermediate_death_rate, intermediate_differentiation_rate, intermediate_dedifferentiation_rate,
differentiated_cell_death_rate, differentiated_cell_dedifferentiation_rate)

# Differentiation vector, ordered from most to least differentiated state (death excluded)
differentiation_rates = [0, intermediate_differentiation_rate, stem_cell_differentiation_rate]

# De-differentiation vector, ordered from most to least differentiated state (death excluded)
dedifferentiation_rates = [differentiated_cell_dedifferentiation_rate, intermediate_dedifferentiation_rate, 0]

# Death vector, ordered from most to least differentiated state
death_rates = [differentiated_cell_death_rate, intermediate_death_rate, stem_cell_death_rate]

# Define transition matrix
# transition_matrix = [1 0 0; stem_cell_death_rate stem_cell_stay_rate stem_cell_differentation_rate; differentiated_cell_death_rate differentiated_cell_dedifferentiation_rate differentiated_cell_stay_rate]

# States: n+3: stem cell; n+2: intermediate state (n); ...; 
# 2: intermediate state (1); 1: differentiated cell; 0: no cell
# where n is the number of intermediate states (an integer)

# Number of simulations to run
number_of_simulations = 5000

# Number of time steps
number_of_time_steps = 751

# Initialise an array to store the chains value over time
state_vector = zeros(number_of_time_steps,number_of_simulations)

# Initialise vector to store extinction times
extinction_times = fill(NaN,number_of_simulations)

# Define initial condition
initial_condition = number_of_intermediate_states + 2 # chain starts with a stem cell

# Plot trajectories
p1 = plot()

# Debugging why things aren't happening as expected
# number_of_deaths_from_S = 0
# number_of_deaths_from_D = 0

for j = 1:number_of_simulations

    global previous_state = initial_condition
    global current_state = initial_condition

    state_vector[1,j] = initial_condition

    for i = 2:number_of_time_steps
        # Run a Markov chain simulation
        u1 = rand(Uniform(0,1))
        if u1 < death_rates[previous_state]
            # Death
            global current_state = 0
            extinction_times[j] = i-1
            break
        elseif u1 > death_rates[previous_state] && u1 < (death_rates[previous_state] + differentiation_rates[previous_state])
            # Differentiation
            global current_state = previous_state - 1
        elseif u1 > (death_rates[previous_state] + differentiation_rates[previous_state]) && 
            u1 < (death_rates[previous_state] + differentiation_rates[previous_state] + dedifferentiation_rates[previous_state])
            # De-differentiation
            global current_state = previous_state + 1
         #else # O/w chain at zero already, will stay there forever
         #    current_state = 0
         #    break
        end
        state_vector[i,j] = current_state
        global previous_state = current_state
    end

    if j <= 10
        plot!(p1,0:number_of_time_steps-1,state_vector[:,j],lw=2,alpha=0.5)
    end

end

# Plot extinction times
# p2 = histogram(extinction_times,color=:turquoise3,label="Simulations",dpi=300,xlimits=(0,number_of_time_steps))
# vline!([expected_extinction_time],lc=:black,lw=2,ls=:dash,label="Mean")
# plot!(fg_legend=:transparent, bg_legend=:transparent)
# xlabel!("Extinction times")
# ylabel!("Frequency")
