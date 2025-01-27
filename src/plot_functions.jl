
@recipe function f(f::BBOBFunction; nlevels = 15, zoom=1)

    markersize := 3
    titlefontsize := 10

    mx = 5.0
    
    if zoom > 1
        x = LinRange(f.x_opt[1]-5/zoom, f.x_opt[1]+5/zoom, 1000)
        y = LinRange(f.x_opt[2]-5/zoom, f.x_opt[2]+5/zoom, 1000)
    else
        x = LinRange(-mx, mx, 1000)
        y = LinRange(-mx, mx, 1000)
    end
    xlim = (minimum(x), maximum(x))
    ylim = (minimum(y), maximum(y))

    z = [f([x,y]) - f.f_opt for y in y, x in x]
    levels = 10 .^(LinRange(-6, log10(maximum(z)), nlevels))

    @series begin
        seriestype := :heatmap
        
        xlabel := "x"
        ylabel := "y"
        cmap := :coolwarm
        alpha := 0.25
        aspectratio := 1
        title := f.name
    
        x, y, z
    end

    @series begin
        seriestype := :contour
        levels := levels
        x, y, z
    end

    @series begin
        seriestype := :scatter
        x = [f.x_opt[1]]
        y = [f.x_opt[2]]
        c := "red"
        label := false
        xlim := xlim
        ylim := ylim

        x,y
    end

end