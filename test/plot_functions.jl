using Gadfly, BBOBFunctions


f = BBOBFunctions.F6

r = 5

plot(
    layer(x=[f.x_opt[1]],y=[f.x_opt[2]],Geom.point,Theme(default_color=colorant"red")),
    layer(
        z=(x,y)->f([x,y])-f.f_opt, x=linspace(-r,r,800), y=linspace(-r,r,800), 
#        Geom.contour(levels=30),
        Geom.contour(levels=logspace(-6,10,40)),
    ),
    Coord.cartesian(xmin=-r,xmax=r,ymin=-r,ymax=r),
)

##
x = collect(linspace(-r,r,1000))

plot(
    layer(x=x,y=[f([xi,0,0]) for xi in x]-f.f_opt,Geom.line),
    layer(x=[f.x_opt[1]],y=[f([f.x_opt[1],0,0])]-f.f_opt,Geom.point),
    layer(x=x,y=[f([0,xi,0]) for xi in x]-f.f_opt,Geom.line,Theme(default_color=colorant"orange")),
    layer(x=[f.x_opt[2]],y=[f([0,f.x_opt[2],0])]-f.f_opt,Geom.point,Theme(default_color=colorant"orange")),
    layer(x=x,y=[f([xi,xi,xi]) for xi in x]-f.f_opt,Geom.line,Theme(default_color=colorant"red")),
    Coord.cartesian(xmin=-r,xmax=r),
    Scale.y_log10
)

#