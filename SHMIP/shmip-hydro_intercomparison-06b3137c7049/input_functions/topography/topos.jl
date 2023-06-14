# Reference surface and bed topography functions, tested with Julia 0.4.6.

## sqrt
"""
Ice sheet margin-like topography.
"""
sqrt_surface(x,y) = 6*( sqrt(x+5e3) - sqrt(5e3) ) + 1
sqrt_bed(x,y) = 0

## Valley
"""
Valley glacier like geometry.
"""
valley_surface(x,y) = 100(x+200)^(1/4) + 1/60*x - 2e10^(1/4) + 1
valley_bed(x,y, para) = f(x,para) + g(y)*h(x,para)

"Parameter to mimic Bench Glacier"
const para_bench = 0.05

# with helper functions
"Thalweg longitudinal profile"
f(x,para) = (valley_surface(6e3,0) - para*6e3)/6e3^2 * x^2 + para*x
"Cross-sectional profile"
g(y) = 0.5e-6 * abs(y)^3
"Correction to cross-sectional profile to keep outline constant"
h(x, para) = (-4.5*x/6e3 + 5) * (valley_surface(x,0)-f(x, para))/(valley_surface(x,0)-f(x, para_bench)+eps())

"""Outline of valley glacier (branch with y>=0"""
valley_outline(x) = ginv( (valley_surface(x,0)-f(x,para_bench))/(h(x,para_bench)+eps()) )
# with inverse of g:
ginv(x) = (x/0.5e-6).^(1/3)

### Suite E
# this is the only suite which varies topography
paras = [0.05, 0, -0.1, -0.5, -0.7]
valley_bed_E(x,y,run) = valley_bed(x,y,paras[run])
valley_surface_E(x,y,run) = valley_surface(x,y)

## Flow-line topographies

sqrt_surface_flowline(x) = sqrt_surface(x,0)
sqrt_bed_flowline(x) = sqrt_bed(x,0)
sqrt_width_flowline(x) = 20e3

valley_surface_flowline(x) = valley_surface(x,0)
valley_bed_flowline(x, para) = f(x,para)
valley_width_flowline(x) = 2*valley_outline(x)

##########
# Plotting
##########
# set to `true` to enable plotting.  Requires the package PyPlot,
# install with `Pkg.add("PyPlot").
if false
    eval(:(using PyPlot))

    # sqrt
    ######
    x = 0:500:100e3
    y = 0:500:20e3
    sx = [sqrt_surface(xx,0) for xx in x]
    bx = [sqrt_bed(xx,0) for xx in x]
    figure(figsize=(8, 4))
    plot(x/1e3, sx)
    hold(true)
    plot(x/1e3, bx)
    xlabel("x (km)")
    ylabel("z (m)")
    title("sqrt topo")

    savefig("sqrt.png")

    # # 3D doesn't work so well
    # s = [sqrt_surface(xx,yy) for xx in x, yy in y]
    # figure()
    # plot_surface(x'/1e3,y''/1e3,s')
    # xlabel("x (km)")
    # ylabel("y (km)")
    # zlabel("z (m)")
    # title("sqrt topo")

    # valley
    ########
    para = para_bench # pick one
    x = 0:20:6e3
    y = -550:20:550
    sx = [valley_surface(xx,0) for xx in x]
    bx = [valley_bed(xx,0,para) for xx in x]

    figure(figsize=(8, 8.5))
    subplot(2,1,1)
    plot(x, sx)
    hold(true)
    plot(x, bx)
    xlabel("x (m)")
    ylabel("z (m)")
    title("valley topo")

    outl = Float64[valley_outline(xx) for xx in x]
    s = Float64[abs(yy)>valley_outline(xx) ? NaN: valley_surface(xx,yy) for xx in x, yy in y]
    b = Float64[abs(yy)>valley_outline(xx) ? NaN: valley_bed(xx,yy,para) for xx in x, yy in y]

    subplot(2,1,2)
    hold(true)
    plot(x, outl, "k", lw=2)
    plot(x, -outl, "k", lw=2)

     # contour is a bit odd, thus all the '' and .+
    im = contour(x'.+0*y'',0*x' .+y'',s', 0:100:1000, colors="b")
    im = contour(x'.+0*y'',0*x' .+y'',b', 0:100:1000, colors="g")

    # Should give outline too:
    # contour(x'.+0*y'',0*x' .+y'',b'-s',[0])
    title("Contours and outline")
    xlabel("x (m)")
    ylabel("y (m)")

    savefig("valley.png")


    # # 3D doesn't work so well
    # figure()
    # hold(true)
    # plot_surface(x',y'',b')
    # plot_surface(x',y'',s')
    # xlabel("x (m)")
    # ylabel("y (m)")
    # zlabel("z (m)")
    # title("valley topo")
end
