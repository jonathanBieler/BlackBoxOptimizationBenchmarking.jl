using Gadfly, BBOBFunctions

##

function plot_fun(f)
    r = 5
    
    p1 = plot(
        layer(x=[f.x_opt[1]],y=[f.x_opt[2]],Geom.point,Theme(default_color=colorant"red")),
        layer(
            z=(x,y)->f([x,y])-f.f_opt, x=linspace(-r,r,800), y=linspace(-r,r,800), 
    #        Geom.contour(levels=30),
            Geom.contour(levels=logspace(-6,10,40)),
        ),
        Coord.cartesian(xmin=-r,xmax=r,ymin=-r,ymax=r),
    )
    
    # line plot
    x = collect(linspace(-r,r,1000))
    
    p2 = plot(
        layer(x=x,y=[ f([xi,f.x_opt[2],f.x_opt[3]]) for xi in x]-f.f_opt, Geom.line),
        layer(x=x,y=[ f([f.x_opt[1],xi,f.x_opt[3]]) for xi in x]-f.f_opt, Geom.line,Theme(default_color=colorant"orange")),
        layer(x=x,y=[ f([f.x_opt[1],f.x_opt[2],xi]) for xi in x]-f.f_opt, Geom.line,Theme(default_color=colorant"red")),
#        layer(x=[f.x_opt[1]],y=[f([f.x_opt[1],0,0])]-f.f_opt,Geom.point),

#        layer(x=[f.x_opt[2]],y=[f([0,f.x_opt[2],0])]-f.f_opt,Geom.point,Theme(default_color=colorant"orange")),
#        layer(x=x,y=[f([xi,xi,xi]) for xi in x]-f.f_opt,Geom.line,Theme(default_color=colorant"red")),
        Coord.cartesian(xmin=-r,xmax=r),
#        Scale.y_log10
    )
    p1,p2
end

include("/Users/jbieler/.julia/v0.5/BBOBFunctions/src/BBOBFunctions.jl")

##

f = BBOBFunctions.F11
p1,p2 = plot_fun(f)

p1


##