using Plots

mutable struct DoublePendule
    pos_x1::Float64
    pos_y1::Float64
    pos_x2::Float64
    pos_y2::Float64
    theta_1::Float64 # Angle 1 au temps t
    theta_2::Float64 # Angle 2 au temps t
    omega_1::Float64 # vitesse angulaire de m1 au temps t
    omega_2::Float64 # vitesse angulaire de m2 au temps t
    m1::Float64 # masse 1
    m2::Float64 # masse 2
    l1::Float64 # longueur 1
    l2::Float64 # longueur 2
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

function eulerstep!(system::DoublePendule, timestep::Float64)::Nothing
    old_omega_1 = copy(system.omega_1)
    old_omega_2 = copy(system.omega_2)

    old_theta_1 = copy(system.theta_1)
    old_theta_2 = copy(system.theta_2)

    acceleration1 = (-system.g*(2*system.m1+system.m2)*sin(system.theta_1)-system.m2*system.g*sin(system.theta_1-2*system.theta_2)-2*sin(system.theta_1-system.theta_2)*system.m2*(system.omega_2^2*system.l2+system.omega_1^2*system.l1*cos(system.theta_1-system.theta_2)))/(system.l1*(2*system.m1+system.m2-system.m2*cos(2*system.theta_1-2*system.theta_2)))
    acceleration2 = (2*sin(system.theta_1-system.theta_2)*(system.omega_1^2*system.l1*(system.m1+system.m2)+system.g*(system.m1+system.m2)*cos(system.theta_1)+system.omega_2^2*system.l2*system.m2*cos(system.theta_1-system.theta_2)))/(system.l2*(2*system.m1+system.m2-system.m2*cos(2*system.theta_1-2*system.theta_2)))

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

function simulate(theta_1::Float64,theta_2::Float64, m1::Float64, m2::Float64, l1::Float64, l2::Float64, g::Float64, number_of_steps::Int, timestep::Float64)
    pos_x1 = l1 * sin(theta_1)
    pos_y1 = -l1 * cos(theta_1)
    pos_x2 = l1 * sin(theta_1) + l2 * sin(theta_2)
    pos_y2 = -l1 * cos(theta_1) - l2 * cos(theta_2)
    system = DoublePendule(pos_x1,pos_y1,pos_x2,pos_y2,theta_1, theta_2, 0, 0, m1, m2, l1, l2, g, [zeros(Float64,number_of_steps) for _ in 1:8]...)

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

function plotSystem(system::DoublePendule, time_iteration::Int)
    Plots.plot([0,0],[0,-2],color=:black, lw=2)
    Plots.scatter!([system.position_x1[time_iteration]],[system.position_y1[time_iteration]],xlim=(-2,2), ylim=(-2,2))
    Plots.scatter!([system.position_x2[time_iteration]],[system.position_y2[time_iteration]])
    Plots.plot!([0,system.position_x1[time_iteration]], [0, system.position_y1[time_iteration]])
    Plots.plot!([system.position_x1[time_iteration],system.position_x2[time_iteration]], [system.position_y1[time_iteration],system.position_y2[time_iteration]])
end

number_of_steps = 100000
timestep = 0.001

system = simulate(3*pi/2, 3*pi/2, 1.0, 0.8, 1.0, 0.8, 9.81, number_of_steps, timestep)

@gif for i in 1:100:number_of_steps
    plotSystem(system, i)
end 