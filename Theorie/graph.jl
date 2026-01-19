using LaTeXStrings
using Plots

"""
    show_pendule_graph(pos1, pos2)

Draw a schema of the double pendulum

### Arguments
- `pos1` : position of m1
- `pos2` : position of m2
"""
function show_pendule_graph(pos1, pos2)
    size = 20
    Plots.plot(xlims=(-2, size), ylims=(-size,2), legend = false, axis=false, grid=false)
    Plots.quiver!([-2,0],[0,-size], quiver=([2+size,0],[0,2+size]), arrow=true, lw=2, color=:black)
    Plots.annotate!(size, -1, text(L"x",10)) 
    Plots.annotate!(1, 2, text(L"y",10)) 
    Plots.quiver!([0,pos1[1]],[0,pos1[2]], quiver=([pos1[1],pos2[1]-pos1[1]],[pos1[2],pos2[2]-pos1[2]]), arrow=true, lw=2, color=:green)
    Plots.scatter!([pos1[1]],[pos1[2]], markersize=20,series_annotations=L"m_1", color=:lightblue)
    Plots.scatter!([pos2[1]],[pos2[2]], markersize=20,series_annotations=L"m_2", color=:lightblue)
    Plots.annotate!(pos1[1]/2-0.5, pos1[2]/2-0.5, text(L"l_1",10,color=:red)) 
    Plots.annotate!(pos1[1]+(pos2[1]-pos1[1])/2-0.5, pos1[2]+(pos2[2]-pos1[2])/2-0.5, text(L"l_2",10,color=:red)) 
    Plots.quiver!([0,0],[0,0], quiver=([pos1[1],pos2[1]],[pos1[2],pos2[2]]), arrow=true, lw=2, color=:black)
    Plots.annotate!(pos1[1]/2+0.5, pos1[2]/2+0.5, text(L"\overrightarrow{r_1}",10,color=:black)) 
    Plots.annotate!(pos2[1]/2+0.5, pos2[2]/2+0.5, text(L"\overrightarrow{r_2}",10,color=:black))

    Plots.annotate!(pos1[1]/6, pos1[2]/3,text(L"\theta_1", 12, color=:orange))
    Plots.plot!([pos1[1], pos1[1]], [pos1[2], -size], linestyle=:dot, color=:red)
    Plots.annotate!(pos1[1]+1, pos1[2]-2,text(L"\theta_2", 12, color=:orange))
end

show_pendule_graph((4,-10),(12,-15))