using Plots
using ProgressMeter

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
    f1::Float64 # frottement pivot 1
    f2::Float64 # frottement pivot 2
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

function simulate(theta_1::Float64,theta_2::Float64, m1::Float64, m2::Float64, l1::Float64, l2::Float64, f1::Float64, f2::Float64, g::Float64, number_of_steps::Int, timestep::Float64)
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

function plotSystem(system::DoublePendule, time_iteration::Int)
    limit = 0.2
    Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=1,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
    Plots.plot!([-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
    Plots.plot!([0,system.position_x1[time_iteration]], [0, system.position_y1[time_iteration]],color=:black)
    Plots.plot!([system.position_x1[time_iteration],system.position_x2[time_iteration]], [system.position_y1[time_iteration],system.position_y2[time_iteration]],color=:black)
    Plots.plot!([system.position_x2[1:time_iteration]],[system.position_y2[1:time_iteration]],color=:cyan)
    Plots.scatter!([system.position_x1[time_iteration]],[system.position_y1[time_iteration]],color=:red)
    Plots.scatter!([system.position_x2[time_iteration]],[system.position_y2[time_iteration]],color=:red)
end

function calc_epot(system::DoublePendule, time_iteration::Int)
    epot = system.m1*system.g*system.position_y1[time_iteration] + system.m2*system.g*system.position_y2[time_iteration]
    return epot
end

function calc_ecin(system::DoublePendule, time_iteration::Int)
    v1 = system.l1 * system.angle_velocity_1[time_iteration]
    v2 = sqrt((system.l1*system.angle_velocity_1[time_iteration])^2 + (system.l2*system.angle_velocity_2[time_iteration])^2 + 2 * system.l1 * system.l2 * system.angle_velocity_1[time_iteration] * system.angle_velocity_2[time_iteration] * cos(system.angles_1[time_iteration] - system.angles_2[time_iteration])) 
    ecin = 1/2 * system.m1 * v1^2 + 1/2 * system.m2 * v2^2
    return ecin
end

function calc_etot(system::DoublePendule, time_iteration::Int)
    return calc_epot(system, time_iteration) + calc_ecin(system, time_iteration)
end

function main()
    nb_secondes = 20
    precision = 10000
    number_of_steps = nb_secondes * precision
    timestep = 1 / (precision)
    FPS = 30
    plot_step = round(Int, 1 / (timestep*FPS))

    steps = collect(1:plot_step:number_of_steps)
    pb = Progress(length(steps))

    angle1 = Float64(3.19559) # rad - sens antihoraire
    angle2 = Float64(3.23186) # rad - sens antihoraire
    m1 = 0.041 # [kg]
    m2 = 0.020 # [kg] - 0.006
    l1 = 0.09174 # [m]
    l2 = 0.06933 # [m]
    f1 = 1.0
    f2 = 1.0
    g = 9.81

    filename = "double_pendule.mp4"

    system = simulate(angle1, angle2, m1, m2, l1, l2,f1, f2, g, number_of_steps, timestep)

    animation = @animate for i in steps
        plotSystem(system, i)
        next!(pb)
    end

    mp4(animation, filename, fps = FPS)

    display(Plots.plot(1:number_of_steps,[calc_ecin(system,i) for i in 1:number_of_steps], title="Energie cinétique", grid = false, legend = false,ylims=(0,0.15), xlabel="temps [s]", ylabel="Energie [J]"))
    display(Plots.plot(1:number_of_steps,[calc_epot(system,i) for i in 1:number_of_steps], title="Energie potentiel", grid = false, legend = false, label = "epot", ylims=(-0.1,0.1), xlabel="temps [s]", ylabel="Energie [J]"))
    display(Plots.plot!(1:number_of_steps,[calc_ecin(system,i) for i in 1:number_of_steps], title="Energie cinétique + potentiel", legend = true, label = "ecin", ylims=(-0.1,0.2), xlabel="temps [s]", ylabel="Energie [J]"))
    display(Plots.plot(1:number_of_steps,[calc_etot(system,i) for i in 1:number_of_steps], title="Energie totale", grid = false, legend = false, ylims=(-0.1,0.1), xlabel="temps [s]", ylabel="Energie [J]"))
end

main()