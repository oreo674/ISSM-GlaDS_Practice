# Testing and plotting the synthetic valley topography.
# See my notebook 4.2.2016
#
# Mimics Bench-glacier
module A
using PyPlot
meshgrid(x,y) = (repmat(x',length(y),1),repmat(y,1,length(x)))

# domain length
xend = 6e3

# surf para
beta = 1/4
s1, s2, sx0, s3 = 100, 100/xend, -200, 1

# bed para
g1 = .5e-6
alpha = 3
f2_example, f20 = -0.7, 300/xend

# surface with ansatz:
#  surf(x,y) = s(x)
surf(x,y) = s(x) .+ 0*y
s(x) = s1*(x-sx0).^beta + s2*x - s1*(-sx0).^beta + s3
s0 = s # default surface

# bed with ansatz:
#  bed(x,y) = f(x) + g(y)*h(x)
bed(x,y,f2=f20) = g(y).*h(x,f2) .+ f(x,f2)
bed0(x,y) = g(y).*h0(x) .+ f0(x) # default bed
g(y) = g1*abs(y).^alpha          # y-dir function
ginv(x) = (x/g1).^(1/alpha)
f(x, f2) = f1(f2)*x.^2 + f2*x # center-line bed, has one free parameter f2

f1(f2) = (s(xend)-f2*xend)/xend^2 # need f(xend) = s(xend)
f0(x) = f(x, f20)
h0(x) = -4.5.*(x./xend) .+ 5

h(x, f2) = h0(x) .* (s(x)-f(x, f2))./(s0(x)-f0(x))

outline0(x) = ginv((s0(x)-f0(x))./h0(x))
outline(x,f2) = ginv((s(x)-f(x,f2))./h(x,f2)) # outline in general

x = 0:10:xend

close("all")
figure()
subplot(3,1,1)
plot(x, s(x))
plot(x, bed0(x,0))
plot(x, bed(x,0, f2_example))

subplot(3,1,2)
plot(x, s(x)-bed0(x,0))
plot(x, s(x)-bed(x,0,f2_example))

subplot(3,1,3)
plot(x, outline0(x))
plot(x, -outline0(x))
plot(x, outline(x,f2_example))
plot(x, -outline(x,f2_example))
axis("equal")

# 3D
y = (-600:10:600)'
X,Y = meshgrid(x,y')

# figure()
# zb0 = bed0(x,y)
# zb = bed(x,y)
# ss = surf(x,y)
# # plot_surface(X,Y,ss')
# # plot_surface(X,Y,zb')
# contour(X,Y,zb0')
# #contour(X,Y,ss')
# plot(x, outline(x), "b", lw=3)
# plot(x, -outline(x), "b", lw=3)
#axis("equal")

nothing
end
