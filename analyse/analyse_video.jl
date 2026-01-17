using VideoIO
using Images
using Statistics
using Plots
using ImageView
using CSV, DataFrames

video_path = "./analyse/video/First_Video_2s.mp4"
filename = "./analyse/double_pendule.mp4"
nframes = 200

l1 = 0.09174 # [m]
l2 = 0.06933 # [m]


# sqrt((492-934)^2+(1872-1884)^2) =~ 443
px_per_meter = 443/0.09174
l1_px = l1 * px_per_meter
l2_px = l2 * px_per_meter

pivot_x = 1884
pivot_y = 934

FPS = 100

function detect_orange_points(frame)
    frame_rgb = RGB.(frame)
    
    r = red.(frame_rgb)
    g = green.(frame_rgb)
    b = blue.(frame_rgb)
    
    # ~rgb(163, 78, 39) => rgb(0.64, 0.31, 0.15) 
    # hsv
    mask = (r .> 0.50) .& (r .< 0.80) .& 
           (g .> 0.15) .& (g .< 0.45) .& 
           (b .> 0.0) .& (b .< 0.30)
    
    return mask
end

function find_mass_center(mask)
    labels = label_components(mask)
    centers = []
    
    for i in 1:maximum(labels)
        coords = findall(labels .== i)
        
        if length(coords) > 35
            y_mean = mean(c[1] for c in coords)
            x_mean = mean(c[2] for c in coords)
            push!(centers, (x_mean, y_mean))
        end
    end
    
    return centers
end

function identify_mass(pos1, pos2, pivot, l1_px)
    d1 = sqrt((pos1[1] - pivot[1])^2 + (pos1[2] - pivot[2])^2)
    d2 = sqrt((pos2[1] - pivot[1])^2 + (pos2[2] - pivot[2])^2)

    m1 = NaN
    m2 = NaN

    if abs(d1 - l1_px) < abs(d2 - l1_px)
        m1 = pos1
        m2 = pos2
    else
        m1 = pos2
        m2 = pos1
    end

    return m1, m2
end

function plotSystem(angle1_list, angle2_list, l1::Float64, l2::Float64, time_iteration::Int, save::Bool)
    x1 = l1 * sin(angle1_list[time_iteration])
    y1 = -l1 * cos(angle1_list[time_iteration])

    x2 = x1 + l2 * sin(angle2_list[time_iteration])
    y2 = y1 - l2 * cos(angle2_list[time_iteration])

    if !save
        limit = 0.2
        Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=1,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
        Plots.plot!([-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
        Plots.plot!([0,x1], [0, y1],color=:black)
        Plots.plot!([x1,x2], [y1,y2],color=:black)
        Plots.scatter!([x1],[y1],color=:red)
        Plots.scatter!([x2],[y2],color=:red)
        return Nothing
    else
        limit = 0.2
        p = Plots.plot([0,0],[0,-limit],color=:black, lw=2, aspect_ratio=1,xlim=(-limit,limit), ylim=(-limit,limit), grid=false, legend=false, axis=false, dpi = 300)
        Plots.plot!(p, [-limit/2,limit/2],[-limit,-limit],color=:black, lw=4)
        Plots.plot!([0,x1], [0, y1],color=:black)
        Plots.plot!([x1,x2], [y1,y2],color=:black)
        Plots.scatter!([x1],[y1],color=:red)
        Plots.scatter!([x2],[y2],color=:red)
        return p
    end
end

function main()
    video = VideoIO.openvideo(video_path)
    angle1_list = []
    angle2_list = []

    for i in 1:nframes
        frame = read(video)
        mask = detect_orange_points(frame)

        centers = find_mass_center(mask)

        pos1, pos2 = centers
        pivot = (pivot_x,pivot_y)
        m1, m2 = identify_mass(pos1, pos2, pivot,l1_px)

        dx1 = m1[1] - pivot_x
        dy1 = m1[2] - pivot_y

        ux1 = 0
        uy1 = 1

        angle1 = acos((dx1*ux1 + dy1*uy1) / (sqrt(dx1^2+dy1^2)*sqrt(ux1^2+uy1^2)))


        if dx1 < 0
            angle1 = 2*pi - angle1
        end



        dx2 = m2[1] - m1[1]
        dy2 = m2[2] - m1[2]

        ux2 = 0
        uy2 = 1

        angle2 = acos((dx2*ux2 + dy2*uy2) / (sqrt(dx2^2+dy2^2)*sqrt(ux2^2+uy2^2)))

        if dx2 < 0
            angle2 = 2*pi - angle2
        end

        push!(angle1_list, angle1)
        push!(angle2_list, angle2)
    end

    df = DataFrame(frame=1:nframes, angle1=angle1_list, angle2=angle2_list)

    CSV.write("./analyse/angles.csv", df)

    animation = @animate for i in 1:1:nframes
        plotSystem(angle1_list, angle2_list, l1, l2, i, false)
    end

    mp4(animation, filename, fps = FPS)
end

main()