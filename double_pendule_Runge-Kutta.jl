using Plots
using ProgressMeter
using CSV
using DataFrames
using Statistics
using Optim

#region System

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

#endregion

#region Simulation

"""
    computeacceleration1(system::DoublePendule,omega_1::Float64,
                         omega_2::Float64,theta_1::Float64,
                         theta_2::Float64) -> Float64

Compute the angular acceleration of the first pendulum.

### Arguments
- `system::DoublePendule` : double pendulum system containing physical parameters
- `omega_1::Float64` : angular velocity of the first pendulum [rad/s]
- `omega_2::Float64` : angular velocity of the second pendulum [rad/s]
- `theta_1::Float64` : angle of the first pendulum [rad]
- `theta_2::Float64` : angle of the second pendulum [rad]

### Returns
- `Float64` : angular acceleration of the first pendulum [rad/s^2]
"""
function computeacceleration1(system::DoublePendule, omega_1::Float64,
                              omega_2::Float64, theta_1::Float64, 
                              theta_2::Float64)

    return (-system.g*(2*system.m1+system.m2)*sin(theta_1)-system.m2*system.g*sin(theta_1-2*theta_2)-2*sin(theta_1-theta_2)*system.m2*(omega_2^2*system.l2+omega_1^2*system.l1*cos(theta_1-theta_2)))/(system.l1*(2*system.m1+system.m2-system.m2*cos(2*theta_1-2*theta_2)))
end

"""
    computeacceleration2(system::DoublePendule,omega_1::Float64,
                         omega_2::Float64,theta_1::Float64,
                         theta_2::Float64) -> Float64

Compute the angular acceleration of the second pendulum.

### Arguments
- `system::DoublePendule` : double pendulum system containing physical parameters
- `omega_1::Float64` : angular velocity of the first pendulum [rad/s]
- `omega_2::Float64` : angular velocity of the second pendulum [rad/s]
- `theta_1::Float64` : angle of the first pendulum [rad]
- `theta_2::Float64` : angle of the second pendulum [rad]

### Returns
- `Float64` : angular acceleration of the first pendulum [rad/s^2]
"""
function computeacceleration2(system::DoublePendule, omega_1::Float64, 
                              omega_2::Float64, theta_1::Float64, 
                              theta_2::Float64)

    return (2*sin(theta_1-theta_2)*(omega_1^2*system.l1*(system.m1+system.m2)+system.g*(system.m1+system.m2)*cos(theta_1)+omega_2^2*system.l2*system.m2*cos(theta_1-theta_2)))/(system.l2*(2*system.m1+system.m2-system.m2*cos(2*theta_1-2*theta_2)))
end

"""
    rungekuttastep!(system::DoublePendule, timestep::Float64) -> Nothing

Advance the state of a double pendulum by one time step using a
fourth-order Runge–Kutta (RK4) integration scheme.

### Arguments
- `system::DoublePendule` : double pendulum system to be advanced in time
- `timestep::Float64` : integration time step Δt [s]

### Sources
- RK4: https://www.youtube.com/watch?v=HOWJp8NV5xU
- Physical equations: 
    - https://www.youtube.com/watch?v=SXj1P9Ra5AM 
    - https://www.myphysicslab.com/pendulum/double-pendulum-en.html 
"""
function rungekuttastep!(system::DoublePendule, timestep::Float64)::Nothing
    acceleration1 = computeacceleration1(system, system.omega_1, system.omega_2, system.theta_1, system.theta_2)
    acceleration2 = computeacceleration2(system, system.omega_1, system.omega_2, system.theta_1, system.theta_2)

    # K1
    k1_omega_1 = acceleration1 - system.f1 * system.omega_1
    k1_omega_2 = acceleration2 - system.f2 * system.omega_2
    k1_theta_1 = system.omega_1
    k1_theta_2 = system.omega_2

    # K2
    omega_1_k2 = system.omega_1 + timestep/2 * k1_omega_1
    omega_2_k2 = system.omega_2 + timestep/2 * k1_omega_2
    theta_1_k2 = system.theta_1 + timestep/2 * k1_theta_1
    theta_2_k2 = system.theta_2 + timestep/2 * k1_theta_2

    k2_omega_1 = computeacceleration1(system, omega_1_k2, omega_2_k2, theta_1_k2, theta_2_k2) - system.f1 * omega_1_k2
    k2_omega_2 = computeacceleration2(system, omega_1_k2, omega_2_k2, theta_1_k2, theta_2_k2) - system.f2 * omega_2_k2
    k2_theta_1 = omega_1_k2
    k2_theta_2 = omega_2_k2

    # K3
    omega_1_k3 = system.omega_1 + timestep/2 * k2_omega_1
    omega_2_k3 = system.omega_2 + timestep/2 * k2_omega_2
    theta_1_k3 = system.theta_1 + timestep/2 * k2_theta_1
    theta_2_k3 = system.theta_2 + timestep/2 * k2_theta_2

    k3_omega_1 = computeacceleration1(system, omega_1_k3, omega_2_k3, theta_1_k3, theta_2_k3) - system.f1 * omega_1_k3
    k3_omega_2 = computeacceleration2(system, omega_1_k3, omega_2_k3, theta_1_k3, theta_2_k3) - system.f2 * omega_2_k3
    k3_theta_1 = omega_1_k3
    k3_theta_2 = omega_2_k3

    # K4
    omega_1_k4 = system.omega_1 + timestep * k3_omega_1
    omega_2_k4 = system.omega_2 + timestep * k3_omega_2
    theta_1_k4 = system.theta_1 + timestep * k3_theta_1
    theta_2_k4 = system.theta_2 + timestep * k3_theta_2

    k4_omega_1 = computeacceleration1(system, omega_1_k4, omega_2_k4, theta_1_k4, theta_2_k4) - system.f1 * omega_1_k4
    k4_omega_2 = computeacceleration2(system, omega_1_k4, omega_2_k4, theta_1_k4, theta_2_k4) - system.f2 * omega_2_k4
    k4_theta_1 = omega_1_k4
    k4_theta_2 = omega_2_k4

    # Update values
    system.omega_1 = system.omega_1 + timestep/6 * (k1_omega_1 + 2*k2_omega_1 + 2*k3_omega_1 + k4_omega_1)
    system.omega_2 = system.omega_2 + timestep/6 * (k1_omega_2 + 2*k2_omega_2 + 2*k3_omega_2 + k4_omega_2)
    
    system.theta_1 = system.theta_1 + timestep/6 * (k1_theta_1 + 2*k2_theta_1 + 2*k3_theta_1 + k4_theta_1)
    system.theta_2 = system.theta_2 + timestep/6 * (k1_theta_2 + 2*k2_theta_2 + 2*k3_theta_2 + k4_theta_2)


    system.pos_x1 = system.l1*sin(system.theta_1)
    system.pos_y1 = -system.l1*cos(system.theta_1)
    system.pos_x2 = system.l1*sin(system.theta_1) + system.l2*sin(system.theta_2)
    system.pos_y2 = -system.l1*cos(system.theta_1) - system.l2*cos(system.theta_2)

    return nothing
end

"""
    simulate(theta_1, theta_2, omega_1, omega_2,
             m1, m2, l1, l2, f1, f2, g,
             number_of_steps, timestep) -> DoublePendule

Simulate the time evolution of a double pendulum.

### Arguments
- `theta_1::Float64` : initial angle of the first pendulum [rad]
- `theta_2::Float64` : initial angle of the second pendulum [rad]
- `omega_1::Float64` : initial angular velocity of the first pendulum [rad/s]
- `omega_2::Float64` : initial angular velocity of the second pendulum [rad/s]

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
function simulate(theta_1::Float64,theta_2::Float64,omega_1::Float64, 
                  omega_2::Float64, m1::Float64, m2::Float64, l1::Float64, 
                  l2::Float64, f1::Float64, f2::Float64, g::Float64, 
                  number_of_steps::Int, timestep::Float64)
    
    

    pos_x1 = l1 * sin(theta_1)
    pos_y1 = -l1 * cos(theta_1)

    pos_x2 = l1 * sin(theta_1) + l2 * sin(theta_2)
    pos_y2 = -l1 * cos(theta_1) - l2 * cos(theta_2)

    system = DoublePendule(pos_x1,pos_y1,pos_x2,pos_y2,theta_1, theta_2, omega_1, omega_2, m1, m2, l1, l2, f1, f2, g, [zeros(Float64,number_of_steps) for _ in 1:8]...)

    system.angles_1[1] = system.theta_1
    system.angles_2[1] = system.theta_2

    system.angle_velocity_1[1] = system.omega_1
    system.angle_velocity_2[1] = system.omega_2

    system.position_x1[1] = system.pos_x1
    system.position_y1[1] = system.pos_y1

    system.position_x2[1] = system.pos_x2
    system.position_y2[1] = system.pos_y2

    for t in 2:number_of_steps
        rungekuttastep!(system, timestep)

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

#endregion

#region Plots

"""
    plotSystem(system, time_iteration, precision) -> Plot


This function generates a plot of the double pendulum system,

### Arguments
- `system::DoublePendule` : simulated double pendulum system
- `time_iteration::Int` : index of the time step to display
- `precision::Int` : number of simulation steps per second

### Returns
- `Plot` : a `Plots.jl` figure showing the pendulum configuration
"""
function plotSystem(system::DoublePendule, time_iteration::Int, precision::Int)
    limit = 0.2
    p = Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=:equal,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
    Plots.plot!(p, [-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
    Plots.plot!(p, [0,system.position_x1[time_iteration]], [0, system.position_y1[time_iteration]],color=:black)
    Plots.plot!(p, [system.position_x1[time_iteration],system.position_x2[time_iteration]], [system.position_y1[time_iteration],system.position_y2[time_iteration]],color=:black)
    Plots.plot!(p, [system.position_x1[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],[system.position_y1[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],color=:cyan)
    Plots.plot!(p, [system.position_x2[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],[system.position_y2[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],color=:cyan)

    Plots.scatter!(p, [system.position_x1[time_iteration]],[system.position_y1[time_iteration]],color=:red)
    Plots.scatter!(p, [system.position_x2[time_iteration]],[system.position_y2[time_iteration]],color=:red)
    return p
end

"""
    plotSystemComparaison(
        system,
        time_iteration,
        precision,
        positions_x1_video,
        positions_y1_video,
        positions_x2_video,
        positions_y2_video
    ) -> Plot

Compare a simulated double pendulum with the tracked video data.

### Arguments
- `system::DoublePendule` : simulated double pendulum system
- `time_iteration::Int` : simulation time step index
- `precision::Int` : simulation steps per second
- `positions_x1_video::Vector{Float64}` : x-coordinates of mass 1 from video
- `positions_y1_video::Vector{Float64}` : y-coordinates of mass 1 from video
- `positions_x2_video::Vector{Float64}` : x-coordinates of mass 2 from video
- `positions_y2_video::Vector{Float64}` : y-coordinates of mass 2 from video

### Returns
- `Plot` : a `Plots.jl` figure comparing simulation and video data
"""
function plotSystemComparaison(system::DoublePendule, time_iteration::Int, precision::Int,
                               positions_x1_video::Vector{Float64}, positions_y1_video::Vector{Float64}, 
                               positions_x2_video::Vector{Float64}, positions_y2_video::Vector{Float64})

    current_index = round(Int, 1 + time_iteration/100)

    x1_video = positions_x1_video[current_index]
    y1_video = positions_y1_video[current_index]

    x2_video = positions_x2_video[current_index]
    y2_video = positions_y2_video[current_index]

    limit = 0.2
    p = Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=:equal,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
    Plots.plot!(p, [-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
    Plots.plot!(p, [0,system.position_x1[time_iteration]], [0, system.position_y1[time_iteration]],color=:black)
    Plots.plot!(p, [system.position_x1[time_iteration],system.position_x2[time_iteration]], [system.position_y1[time_iteration],system.position_y2[time_iteration]],color=:black)
    Plots.plot!(p, [system.position_x1[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],[system.position_y1[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],color=:cyan)
    Plots.plot!(p, [system.position_x2[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],[system.position_y2[max(time_iteration-round(Int,0.5*precision),1):time_iteration]],color=:cyan)    
    Plots.scatter!(p, [system.position_x1[time_iteration]],[system.position_y1[time_iteration]],color=:red)
    Plots.scatter!(p, [system.position_x2[time_iteration]],[system.position_y2[time_iteration]],color=:red)

    Plots.plot!(p, [positions_x1_video[max(current_index-25,1):current_index]],[positions_y1_video[max(current_index-25,1):current_index]],color=:yellow)
    Plots.plot!(p, [positions_x2_video[max(current_index-25,1):current_index]],[positions_y2_video[max(current_index-25,1):current_index]],color=:yellow)
    Plots.scatter!(p,[x1_video],[y1_video],color=:blue)
    Plots.scatter!(p,[x2_video],[y2_video],color=:blue)
    return p
end

"""
    plotSystemPositions(
        system,
        precision,
        positions_x1_video,
        positions_y1_video,
        positions_x2_video,
        positions_y2_video
    )

Plot simulated and video-tracked positions.

### Arguments
- `system::DoublePendule` : simulated double pendulum system
- `precision::Int` : number of simulation steps per second
- `positions_x1_video::Vector{Float64}` : x-coordinates of mass 1 from video tracking
- `positions_y1_video::Vector{Float64}` : y-coordinates of mass 1 from video tracking
- `positions_x2_video::Vector{Float64}` : x-coordinates of mass 2 from video tracking
- `positions_y2_video::Vector{Float64}` : y-coordinates of mass 2 from video tracking
"""
function plotSystemPositions(system::DoublePendule, precision::Int, positions_x1_video::Vector{Float64}, 
                             positions_y1_video::Vector{Float64}, positions_x2_video::Vector{Float64}, 
                             positions_y2_video::Vector{Float64})

    nTime = round(Int, precision / 100)

    p = plot(layout = (2, 2))

    Plots.plot!(p[1], 1:length(system.position_x1), [system.position_x1], title="x1", grid=false, legend=false, xlabel="time", ylabel="values")
    Plots.plot!(p[1], 1:length(system.position_x1), [repeat(positions_x1_video, inner = nTime)])

    Plots.plot!(p[2], 1:length(system.position_y1), [system.position_y1], title="y1", grid=false, legend=false, xlabel="time", ylabel="values")
    Plots.plot!(p[2], 1:length(system.position_y1), [repeat(positions_y1_video, inner = nTime)])

    Plots.plot!(p[3], 1:length(system.position_x2), [system.position_x2], title="x2", grid=false, legend=false, xlabel="time", ylabel="values")
    Plots.plot!(p[3], 1:length(system.position_x2), [repeat(positions_x2_video, inner = nTime)])

    Plots.plot!(p[4], 1:length(system.position_y2), [system.position_y2], title="y2", grid=false, legend=false, xlabel="time", ylabel="values")
    Plots.plot!(p[4], 1:length(system.position_y2), [repeat(positions_y2_video, inner = nTime)])
    
    display(p)
end

#endregion

#region Energie

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
    displayEnergie(system::DoublePendule, number_of_steps::Int)

Display the energy evolution of a double pendulum over time using plots.

### Arguments
- `system::DoublePendule` : the simulated double pendulum system
- `number_of_steps::Int` : the number of simulation steps to display
"""
function displayEnergie(system::DoublePendule, number_of_steps::Int)
    display(Plots.plot(1:number_of_steps,[calc_ecin(system,i) for i in 1:number_of_steps], title="Energie cinétique", grid = false, legend = false,ylims=(0,0.15), xlabel="temps [ms]", ylabel="Energie [J]"))
    display(Plots.plot(1:number_of_steps,[calc_epot(system,i) for i in 1:number_of_steps], title="Energie potentiel", grid = false, legend = false, label = "epot", ylims=(-0.1,0.1), xlabel="temps [s]", ylabel="Energie [J]"))
    display(Plots.plot!(1:number_of_steps,[calc_ecin(system,i) for i in 1:number_of_steps], title="Energie cinétique + potentiel", legend = true, label = "ecin", ylims=(-0.1,0.2), xlabel="temps [s]", ylabel="Energie [J]"))
    display(Plots.plot(1:number_of_steps,[calc_etot(system,i) for i in 1:number_of_steps], title="Energie totale", grid = false, legend = false, ylims=(-0.1,0.1), xlabel="temps [ms]", ylabel="Energie [J]"))
end

#endregion

#region Find parameters

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
    error = deg2rad(0.5)

    diff = abs.(angle_video .- angle_sim) .- error

    diff = max.(diff, 0.0)  

    return sqrt(mean(diff.^2)) / std(angle_video)
end
#endregion

#region Main

"""
    main(display_video::Bool, display_energie::Bool, find_parameters::Bool)

Run the full double pendulum simulation and optionally display video, energy plots, and optimize parameters.

### Arguments
- `display_video::Bool` : if true, creates an animation of the double pendulum
- `display_energie::Bool` : if true, plots kinetic, potential, and total energy over time
- `find_parameters::Bool` : if true, performs parameter optimization to fit the simulation to video data
"""
function main(display_video::Bool, display_energie::Bool, find_parameters::Bool)
    nb_secondes = 2 # [s]
    precision = 10000
    number_of_steps = nb_secondes * precision
    timestep = 1 / (precision)
    FPS = 100
    plot_step = round(Int, 1 / (timestep*FPS))

    # Progress bar
    steps = collect(1:plot_step:number_of_steps)
    pb = Progress(length(steps))

    # Parameters
    angle1 = Float64(3.1686318134763938) # [rad]
    angle2 = Float64(3.237188742291009) # [rad]
    omega_1 = 0.0 # [m*s^-1] velocity
    omega_2 = 0.061 # [m*s^-1] velocity
    m1 = 0.021 # [kg]
    m2 = 0.0028 # [kg]
    l1 = 0.09174 # [m]
    l2 = 0.06933 # [m]
    f1 = 3.211948814733922e-5 # 0.15 [s^-1]
    f2 = 6.133546094892238e-6 # 0.15 [s^-1]
    g = 9.81

    m_factor = 6.857892359965417
    m1 = 1.0
    m2 = m1/m_factor

    # Video
    df = CSV.read("./analyse/angles.csv", DataFrame)
    number_of_frames = 20

    angles1_video = Vector{Float64}(df.angle1)
    angles2_video = Vector{Float64}(df.angle2)

    positions_x1_video = [l1 * sin(angle1_video) for angle1_video in angles1_video]
    positions_y1_video = [-l1 * cos(angle1_video) for angle1_video in angles1_video]

    positions_x2_video = [l1 * sin(angles1_video[i]) + l2 * sin(angles2_video[i]) for i in 1:length(angles2_video)]
    positions_y2_video = [-l1 * cos(angles1_video[i]) - l2 * cos(angles2_video[i]) for i in 1:length(angles2_video)]

    # Simulate
    system = simulate(angle1, angle2, omega_1, omega_2, m1, m2, l1, l2,f1, f2, g, number_of_steps, timestep)

    plotSystemPositions(system, precision, positions_x1_video, positions_y1_video, positions_x2_video, positions_y2_video)

    if display_video
        filename = "double_pendule_Runge-Kutta.mp4"

        animation = @animate for i in steps
            # plotSystem(system, i, precision)
            plotSystemComparaison(system, i, precision, positions_x1_video, positions_y1_video, positions_x2_video, positions_y2_video)

            next!(pb)
        end

        mp4(animation, filename, fps = FPS)
    end

    if display_energie
        displayEnergie(system, number_of_steps)
    end

    if find_parameters
        
        p0_init = [angle1, angle2, omega_1, omega_2]
        lower_init = [pi, pi, -0.1, -0.1]
        upper_init = [pi + pi/4, pi + pi/4, 0.1, 0.1]

        function loss_init(p)
            theta_1, theta_2, omega_1, omega_2 = p

            system = simulate(theta_1, theta_2, omega_1, omega_2, 1.0, 1.0/m_factor, l1, l2, f1, f2, g, number_of_steps, timestep)
            
            step = round(Int, 0.01 / timestep)
            angle1_sim = system.angles_1[1:step:end][1:number_of_frames]
            angle2_sim = system.angles_2[1:step:end][1:number_of_frames]
            
            e1 = NRMSE(angles1_video[1:number_of_frames], angle1_sim)
            e2 = NRMSE(angles2_video[1:number_of_frames], angle2_sim)
            
            return (e1 + e2)/2
        end

        result_init = optimize(loss_init, lower_init, upper_init, p0_init, Fminbox(NelderMead()))
        p_init = Optim.minimizer(result_init)

        angle1, angle2, omega_1, omega_2 = p_init
        println("Phase 1 : conditions initiales optimisées -> ", p_init)
        
        p0_global = [f1, f2, m_factor]
        lower_global = [0.0, 0.0, 2.0]
        upper_global = [0.3, 0.3, 10.0]

        function loss_global(p)
            f1_opt, f2_opt, m_factor_opt = p
            system = simulate(angle1, angle2, omega_1, omega_2, 1.0, 1.0/m_factor_opt, l1, l2, f1_opt, f2_opt, g, number_of_steps, timestep)
            
            step = round(Int, 0.01 / timestep)
            angle1_sim = system.angles_1[1:step:end][1:number_of_frames]
            angle2_sim = system.angles_2[1:step:end][1:number_of_frames]
            
            e1 = NRMSE(angles1_video[1:number_of_frames], angle1_sim)
            e2 = NRMSE(angles2_video[1:number_of_frames], angle2_sim)
            
            return (e1 + e2)/2
        end

        result_global = optimize(loss_global, lower_global, upper_global, p0_global, Fminbox(NelderMead()))
        f1, f2, m_factor = Optim.minimizer(result_global)

        println("Phase 2 : paramètres globaux optimisés -> ", f1, " ", f2, " ", m_factor)
        println("Erreur finale = ", Optim.minimum(result_global))

        system = simulate(angle1, angle2, omega_1, omega_2, 1.0, 1.0/m_factor, l1, l2,f1, f2, g, number_of_steps, timestep)

        plotSystemPositions(system, precision, positions_x1_video, positions_y1_video, positions_x2_video, positions_y2_video)
    end
end

# display_video, display_energie, find_parameters
main(false, false, true)

#endregion