using Gadfly

using BlackBoxOptimizationBenchmarking
const BBOB = BlackBoxOptimizationBenchmarking

##

function plot_fun(f)
    r = 5
    
    c, n, nline = Geom.contour(levels=exp10.(range(-6, stop=10, length=40))), 800, 1000
    if any(f .== [BBOB.F16,BBOB.F17,BBOB.F18,BBOB.F19])
        c, n, nline = Geom.contour(levels=exp10.(range(-6, stop=8, length=20))), 500, 1000
    end

    p1 = plot(
        layer(x=[f.x_opt[1]],y=[f.x_opt[2]],Geom.point,Theme(default_color=colorant"red")),
        layer(
            z=(x,y)->f([x,y])-f.f_opt, x=range(-r, stop=r, length=n), y=range(-r, stop=r, length=n), 
            c,
        ),
        Coord.cartesian(xmin=-r,xmax=r,ymin=-r,ymax=r),
        Guide.title(string(f))
    )
    
    # line plot
    x = collect(range(-r, stop=r, length=nline))
    
    p2 = plot(
        layer(x=x,y=[ f([xi,f.x_opt[2],f.x_opt[3]]) for xi in x]-f.f_opt, Geom.line),
        layer(x=x,y=[ f([f.x_opt[1],xi,f.x_opt[3]]) for xi in x]-f.f_opt, Geom.line,Theme(default_color=colorant"orange")),
        layer(x=x,y=[ f([f.x_opt[1],f.x_opt[2],xi]) for xi in x]-f.f_opt, Geom.line,Theme(default_color=colorant"red")),
        Coord.cartesian(xmin=-r,xmax=r),
        Guide.title(string(f)),
        Scale.y_log10
    )
    p1,p2
end

outdir = joinpath(Pkg.dir(),"BBOB","data","plots","functions")

for i in 20:length(enumerate(BBOB.BBOBFunction))
    
    f = getfield(BBOB,Symbol("F$i"))    
    info(f)
    p1,p2 = plot_fun(f)

    draw(PDF(joinpath(outdir,"F$i.pdf"),16cm,16cm),p1)
    draw(SVG(joinpath(outdir,"F$i.svg"),16cm,16cm),p1)
    draw(PDF(joinpath(outdir,"1D","F$i.pdf"),16cm,12cm),p2)
    draw(SVG(joinpath(outdir,"1D","F$i.svg"),16cm,12cm),p2)
end

##