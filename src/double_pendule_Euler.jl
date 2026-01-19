using Plots
using ProgressMeter
using CSV
using DataFrames
using Statistics

"""
    DoublePendule

Class for storing information about double pendulums.

### Fields
- `pos_x1::Float64` : x-position of the first mass at time *t* [m]
- `pos_y1::Float64` : y-position of the first mass at time *t* [m]
- `pos_x2::Float64` : x-position of the second mass at time *t* [m]
- `pos_y2::Float64` : y-position of the second mass at time *t* [m]

- `theta_1::Float64` : angle of the first pendulum at time *t* [rad]
- `theta_2::Float64` : angle of the second pendulum at time *t* [rad]
- `omega_1::Float64` : angular velocity of the first pendulum at time *t* [rad/s]
- `omega_2::Float64` : angular velocity of the second pendulum at time *t* [rad/s]

- `m1::Float64` : mass of the first pendulum [kg]
- `m2::Float64` : mass of the second pendulum [kg]
- `l1::Float64` : length of the first rod [m]
- `l2::Float64` : length of the second rod [m]

- `f1::Float64` : coefficient of friction at the first pivot [1/s]
- `f2::Float64` : coefficient of friction at the second pivot [1/s]
- `g::Float64`  : gravitational acceleration [m/s^2]

- `angles_1::Vector{Float64}` : history of angle 1
- `angles_2::Vector{Float64}` : history of angle 2

- `angle_velocity_1::Vector{Float64}` : history of angle 1 velocity
- `angle_velocity_2::Vector{Float64}` : history of angle 2 velocity

- `position_x1::Vector{Float64}` : x-position history of the first mass
- `position_y1::Vector{Float64}` : y-position history of the first mass

- `position_x2::Vector{Float64}` : x-position history of the second mass
- `position_y2::Vector{Float64}` : y-position history of the second mass
"""
mutable struct DoublePendule
    pos_x1::Float64
    pos_y1::Float64
    pos_x2::Float64
    pos_y2::Float64
    theta_1::Float64
    theta_2::Float64
    omega_1::Float64
    omega_2::Float64
    m1::Float64
    m2::Float64
    l1::Float64
    l2::Float64
    f1::Float64
    f2::Float64
    g::Float64

    angles_1::Vector{Float64}
    angles_2::Vector{Float64}

    angle_velocity_1::Vector{Float64}
    angle_velocity_2::Vector{Float64}

    position_x1::Vector{Float64}
    position_y1::Vector{Float64}

    position_x2::Vector{Float64}
    position_y2::Vector{Float64}
end

"""
    eulerstep!(system::DoublePendule, timestep::Float64) -> Nothing

Advance the state of a double pendulum by one time step using a
euler integration scheme.

### Arguments
- `system::DoublePendule` : double pendulum system to be advanced in time
- `timestep::Float64` : integration time step Δt [s]

### Sources
- Physical equations: 
    - https://www.youtube.com/watch?v=SXj1P9Ra5AM 
    - https://www.myphysicslab.com/pendulum/double-pendulum-en.html 
"""
function eulerstep!(system::DoublePendule, timestep::Float64)::Nothing
    old_omega_1 = copy(system.omega_1)
    old_omega_2 = copy(system.omega_2)

    old_theta_1 = copy(system.theta_1)
    old_theta_2 = copy(system.theta_2)

    acceleration1 = (-system.g*(2*system.m1+system.m2)*sin(system.theta_1)-system.m2*system.g*sin(system.theta_1-2*system.theta_2)-2*sin(system.theta_1-system.theta_2)*system.m2*(system.omega_2^2*system.l2+system.omega_1^2*system.l1*cos(system.theta_1-system.theta_2)))/(system.l1*(2*system.m1+system.m2-system.m2*cos(2*system.theta_1-2*system.theta_2)))
    acceleration2 = (2*sin(system.theta_1-system.theta_2)*(system.omega_1^2*system.l1*(system.m1+system.m2)+system.g*(system.m1+system.m2)*cos(system.theta_1)+system.omega_2^2*system.l2*system.m2*cos(system.theta_1-system.theta_2)))/(system.l2*(2*system.m1+system.m2-system.m2*cos(2*system.theta_1-2*system.theta_2)))

    acceleration1 -= system.f1 * system.omega_1
    acceleration2 -= system.f2 * system.omega_2

    system.omega_1 = old_omega_1 + timestep * acceleration1
    system.omega_2 = old_omega_2 + timestep * acceleration2

    system.theta_1 = old_theta_1 + timestep * system.omega_1
    system.theta_2 = old_theta_2 + timestep * system.omega_2

    system.pos_x1 = system.l1*sin(system.theta_1)
    system.pos_y1 = -system.l1*cos(system.theta_1)
    system.pos_x2 = system.l1*sin(system.theta_1) + system.l2*sin(system.theta_2)
    system.pos_y2 = -system.l1*cos(system.theta_1) - system.l2*cos(system.theta_2)

    return nothing
end

"""
    simulate(theta_1, theta_2,
             m1, m2, l1, l2, f1, f2, g,
             number_of_steps, timestep) -> DoublePendule

Simulate the time evolution of a double pendulum.

### Arguments
- `theta_1::Float64` : initial angle of the first pendulum [rad]
- `theta_2::Float64` : initial angle of the second pendulum [rad]

- `m1::Float64` : mass of the first pendulum [kg]
- `m2::Float64` : mass of the second pendulum [kg]
- `l1::Float64` : length of the first rod [m]
- `l2::Float64` : length of the second rod [m]

- `f1::Float64` : coefficient of friction at the first pivot [1/s]
- `f2::Float64` : coefficient of friction at the second pivot [1/s]
- `g::Float64`  : gravitational acceleration [m/s^2]

- `number_of_steps::Int` : number of integration time steps
- `timestep::Float64` : integration time step dt [s]

### Returns
- `DoublePendule` : simulated system containing the full time history
"""
function simulate(theta_1::Float64,theta_2::Float64, 
                  m1::Float64, m2::Float64, l1::Float64, 
                  l2::Float64, f1::Float64, f2::Float64, 
                  g::Float64, number_of_steps::Int, 
                  timestep::Float64)
    
                  pos_x1 = l1 * sin(theta_1)
    pos_y1 = -l1 * cos(theta_1)
    pos_x2 = l1 * sin(theta_1) + l2 * sin(theta_2)
    pos_y2 = -l1 * cos(theta_1) - l2 * cos(theta_2)
    system = DoublePendule(pos_x1,pos_y1,pos_x2,pos_y2,theta_1, theta_2, 0, 0, m1, m2, l1, l2, f1, f2, g, [zeros(Float64,number_of_steps) for _ in 1:8]...)

    system.angles_1[1] = system.theta_1
    system.angles_2[1] = system.theta_2

    system.angle_velocity_1[1] = system.omega_1
    system.angle_velocity_2[1] = system.omega_2

    system.position_x1[1] = system.pos_x1
    system.position_y1[1] = system.pos_y1

    system.position_x2[1] = system.pos_x2
    system.position_y2[1] = system.pos_y2

    for t in 2:number_of_steps
        eulerstep!(system, timestep)

        system.angles_1[t] = system.theta_1
        system.angles_2[t] = system.theta_2

        system.angle_velocity_1[t] = system.omega_1
        system.angle_velocity_2[t] = system.omega_2

        system.position_x1[t] = system.pos_x1
        system.position_y1[t] = system.pos_y1

        system.position_x2[t] = system.pos_x2
        system.position_y2[t] = system.pos_y2
    end

    return system
end

"""
    plotSystem(system, time_iteration, save) -> Plot


This function generates a plot of the double pendulum system,

### Arguments
- `system::DoublePendule` : simulated double pendulum system
- `time_iteration::Int` : index of the time step to display
- `save::Bool` : if true, return the plot

### Returns
- `Plot` : a `Plots.jl` figure showing the pendulum configuration
"""
function plotSystem(system::DoublePendule, time_iteration::Int, save::Bool)
    if !save
        limit = 0.2
        Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=:equal,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
        Plots.plot!([-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
        Plots.plot!([0,system.position_x1[time_iteration]], [0, system.position_y1[time_iteration]],color=:black)
        Plots.plot!([system.position_x1[time_iteration],system.position_x2[time_iteration]], [system.position_y1[time_iteration],system.position_y2[time_iteration]],color=:black)
        Plots.plot!([system.position_x2[max(time_iteration-50,1):time_iteration]],[system.position_y2[max(time_iteration-50,1):time_iteration]],color=:cyan)
        Plots.plot!([system.position_x1[max(time_iteration-50,1):time_iteration]],[system.position_y1[max(time_iteration-50,1):time_iteration]],color=:red)
        Plots.scatter!([system.position_x1[time_iteration]],[system.position_y1[time_iteration]],color=:red)
        Plots.scatter!([system.position_x2[time_iteration]],[system.position_y2[time_iteration]],color=:red)
        return Nothing
    else
        limit = 0.2
        p = Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=:equal,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
        Plots.plot!(p, [-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
        Plots.plot!(p, [0,system.position_x1[time_iteration]], [0, system.position_y1[time_iteration]],color=:black)
        Plots.plot!(p, [system.position_x1[time_iteration],system.position_x2[time_iteration]], [system.position_y1[time_iteration],system.position_y2[time_iteration]],color=:black)
        Plots.plot!(p, [system.position_x2[max(time_iteration-50,1):time_iteration]],[system.position_y2[max(time_iteration-50,1):time_iteration]],color=:cyan)
        Plots.plot!(p, [system.position_x1[max(time_iteration-50,1):time_iteration]],[system.position_y1[max(time_iteration-50,1):time_iteration]],color=:red)
        Plots.scatter!(p, [system.position_x1[time_iteration]],[system.position_y1[time_iteration]],color=:red)
        Plots.scatter!(p, [system.position_x2[time_iteration]],[system.position_y2[time_iteration]],color=:red)
        return p
    end
end

"""
    calc_epot(system::DoublePendule, time_iteration::Int) -> Float64

Compute the potential energy of a double pendulum at a specific time step.

### Arguments
- `system::DoublePendule` : the simulated double pendulum system
- `time_iteration::Int` : the index of the time step at which to compute the energy

### Returns
- `Float64` : total gravitational potential energy [J] at the given time step
"""
function calc_epot(system::DoublePendule, time_iteration::Int)
    epot = system.m1*system.g*system.position_y1[time_iteration] + system.m2*system.g*system.position_y2[time_iteration]
    return epot
end

"""
    calc_ecin(system::DoublePendule, time_iteration::Int) -> Float64

Compute the kinetic energy of a double pendulum at a specific time step.

### Arguments
- `system::DoublePendule` : the simulated double pendulum system
- `time_iteration::Int` : the index of the time step at which to compute the energy

### Returns
- `Float64` : total kinetic energy [J] at the given time step
"""
function calc_ecin(system::DoublePendule, time_iteration::Int)
    v1 = system.l1 * system.angle_velocity_1[time_iteration]
    v2 = sqrt((system.l1*system.angle_velocity_1[time_iteration])^2 + (system.l2*system.angle_velocity_2[time_iteration])^2 + 2 * system.l1 * system.l2 * system.angle_velocity_1[time_iteration] * system.angle_velocity_2[time_iteration] * cos(system.angles_1[time_iteration] - system.angles_2[time_iteration])) 
    ecin = 1/2 * system.m1 * v1^2 + 1/2 * system.m2 * v2^2
    return ecin
end

"""
    calc_etot(system::DoublePendule, time_iteration::Int) -> Float64

Compute the total mechanical energy of a double pendulum at a specific time step.

### Arguments
- `system::DoublePendule` : the simulated double pendulum system
- `time_iteration::Int` : the index of the time step at which to compute the total energy

### Returns
- `Float64` : total mechanical energy [J] at the given time step
"""
function calc_etot(system::DoublePendule, time_iteration::Int)
    return calc_epot(system, time_iteration) + calc_ecin(system, time_iteration)
end

"""
    NRMSE(angle_video::Vector{Float64}, angle_sim::Vector{Float64}) -> Float64

Compute the normalized root mean square error (NRMSE) between simulated angles and video angles.

### Arguments
- `angle_video::Vector{Float64}` : angles from video [rad]
- `angle_sim::Vector{Float64}` : simulated angles [rad]

### Returns
- `Float64` : normalized RMSE
"""
function NRMSE(angle_video::Vector{Float64}, angle_sim::Vector{Float64})
    return sqrt(mean((angle_video .- angle_sim).^2)) / std(angle_video)
end

"""
    find_mass_factor(m_factor::Float64, angle1_video::Vector{Float64}, angle2_video::Vector{Float64},
                     angle1::Float64, angle2::Float64, l1::Float64, l2::Float64, f1::Float64,
                     f2::Float64, g::Float64, timestep::Float64, number_of_steps::Int, number_of_frames::Int)

Compute a score comparing the double pendulum simulation to reference video angles for a given mass factor.

### Arguments
- `m_factor::Float64` : ratio of mass m1 to m2
- `angle1_video::Vector{Float64}` : reference angles of the first pendulum from video [rad]
- `angle2_video::Vector{Float64}` : reference angles of the second pendulum from video [rad]
- `angle1::Float64` : initial angle of the first pendulum [rad]
- `angle2::Float64` : initial angle of the second pendulum [rad]
- `l1::Float64` : length of the first pendulum rod [m]
- `l2::Float64` : length of the second pendulum rod [m]
- `f1::Float64` : coefficient of friction at the first pivot [1/s]
- `f2::Float64` : coefficient of friction at the second pivot [1/s]
- `g::Float64` : gravitational acceleration [m/s^2]
- `timestep::Float64` : simulation time step [s]
- `number_of_steps::Int` : total number of simulation steps
- `number_of_frames::Int` : number of frames to sample for error calculation

### Returns
- `score::Float64` : average NRMSE between simulation and video for both angles
- `nrmse1::Float64` : NRMSE for angle1
- `nrmse2::Float64` : NRMSE for angle2
"""
function find_mass_factor(m_factor::Float64,angle1_video::Vector{Float64},
                          angle2_video::Vector{Float64},angle1::Float64,
                          angle2::Float64,l1::Float64,l2::Float64,f1::Float64,
                          f2::Float64,g::Float64,timestep::Float64, 
                          number_of_steps::Int, number_of_frames::Int)

    m1 = 1.0
    m2 = m1 / m_factor

    system = simulate(angle1,angle2,m1,m2,l1,l2,f1,f2,g,number_of_steps,timestep)
    step = round(Int, 0.01 / timestep)
    
    angle1_system::Vector{Float64} = []
    angle2_system::Vector{Float64} = []

    for i in 1:step:length(system.angles_1)
        push!(angle1_system, system.angles_1[i])
        push!(angle2_system, system.angles_2[i])
    end

    angle1_system = angle1_system[1:number_of_frames]
    angle2_system = angle2_system[1:number_of_frames]

    nrmse1 = NRMSE(angle1_video, angle1_system)
    nrmse2 = NRMSE(angle2_video, angle2_system)

    score = (nrmse1 + nrmse2) / 2

    return score, nrmse1, nrmse2
end

"""
    main(display_video::Bool, display_energie::Bool, find_parameters::Bool)

Run the full double pendulum simulation and optionally display video, energy plots, and optimize parameters.

### Arguments
- `display_video::Bool` : if true, creates an animation of the double pendulum
- `display_energie::Bool` : if true, plots kinetic, potential, and total energy over time
- `find_mass::Bool` : if true, performs parameter optimization to fit the simulation to video data
"""
function main(display_video::Bool, display_energie::Bool, find_mass::Bool)
    nb_secondes = 2 # [s]
    precision = 10000
    number_of_steps = nb_secondes * precision
    timestep = 1 / (precision)
    FPS = 30
    plot_step = round(Int, 1 / (timestep*FPS))

    steps = collect(1:plot_step:number_of_steps)
    pb = Progress(length(steps))

    angle1 = Float64(3.1686318134763938) # rad - 3.1686318134763938 +- 0.0349066 (2°) -> [3.13373, 3.20354]
    angle2 = Float64(3.237188742291009) # rad - 3.237188742291009 +- 0.0349066 (2°) -> [3.20228, 3.2721]
    m1 = 0.021 # [kg]
    m2 = 0.0028 # [kg]
    l1 = 0.09174 # [m]
    l2 = 0.06933 # [m]
    f1 = 0.0 # 0.15 [s^-1]
    f2 = 0.0 # 0.15 [s^-1]
    g = 9.81

    # m1 = 1.0
    # m2 = m1 / 15.40943


    if find_mass
        for angle1 in range(3.13373, 3.20354; length=10)
            for angle2 in range(3.20228, 3.2721; length=10)
        
                df = CSV.read("./analyse/angles.csv", DataFrame)
                number_of_frames = 10

                angles1_video = Vector{Float64}(df.angle1)
                angles2_video = Vector{Float64}(df.angle2)

                angle1_video = angles1_video[1:number_of_frames]
                angle2_video = angles2_video[1:number_of_frames]

                timestep = 1/precision
                # angle1 = angle1_video[1]
                # angle2 = angle2_video[1]

                factors = 10.0:0.001:20.0
                errors = Float64[]
                NRMSE_angles1 = Float64[]
                NRMSE_angles2 = Float64[]

                for m_factor in factors
                    e, NRMSE_angle1, NRMSE_angle2 = find_mass_factor(m_factor,angle1_video,angle2_video,angle1,angle2,l1,l2,f1,f2,g,timestep, number_of_steps, number_of_frames)
                    push!(errors, e)
                    push!(NRMSE_angles1, NRMSE_angle1)
                    push!(NRMSE_angles2, NRMSE_angle2)
                end

                best_masse_idx = argmin(errors)
                println("The best m factor is ", factors[best_masse_idx], " with a NRMSE score of : ", errors[best_masse_idx], " NRMSE score of angle1 is ", NRMSE_angles1[best_masse_idx], " NRMSE score of angle2 is ", NRMSE_angles2[best_masse_idx], " angle1: ", angle1, "angle2: ", angle2)
            
            end
        end
    end

    if display_video
        filename = "src/double_pendule_Euler.mp4"

        system = simulate(angle1, angle2, m1, m2, l1, l2,f1, f2, g, number_of_steps, timestep)

        animation = @animate for i in steps
            plotSystem(system, i, false)
            next!(pb)
        end

        mp4(animation, filename, fps = FPS)
    end

    if display_energie
        display(Plots.plot(1:number_of_steps,[calc_ecin(system,i) for i in 1:number_of_steps], title="Energie cinétique", grid = false, legend = false,ylims=(0,0.15), xlabel="temps [ms]", ylabel="Energie [J]"))
        display(Plots.plot(1:number_of_steps,[calc_epot(system,i) for i in 1:number_of_steps], title="Energie potentiel", grid = false, legend = false, label = "epot", ylims=(-0.1,0.1), xlabel="temps [s]", ylabel="Energie [J]"))
        display(Plots.plot!(1:number_of_steps,[calc_ecin(system,i) for i in 1:number_of_steps], title="Energie cinétique + potentiel", legend = true, label = "ecin", ylims=(-0.1,0.2), xlabel="temps [s]", ylabel="Energie [J]"))
        display(Plots.plot(1:number_of_steps,[calc_etot(system,i) for i in 1:number_of_steps], title="Energie totale", grid = false, legend = false, ylims=(-0.1,0.1), xlabel="temps [ms]", ylabel="Energie [J]"))
    end
end

# display_video, display_energie, find_mass
main(true, false, false)