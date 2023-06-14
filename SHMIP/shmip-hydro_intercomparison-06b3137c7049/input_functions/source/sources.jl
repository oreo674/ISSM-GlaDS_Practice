# Reference implementation of source functions, tested with Julia 0.4.6.
#
# Distributed sources have signature  source_X(x,y,t,run)
# Moulin sources have signature       source_X_m(moulin_id,t,run)
# where X is the suite-letter.
#
# x,y : coordinates in meter
# t : in seconds
# moulin_id : the index of the moulin, starting from 0
# run : run number in a suite

# A #######

"""
Steady distributed source of Suite A, in m/s.
"""
source_A(x,y,t,run) = [7.93e-11, 1.59e-09, 5.79e-09,  2.5e-08,  4.5e-08,  5.79e-07][run]
source_A_m(moulin_id,t,run) = 0

# B #######

const moulins = Matrix{Float64}[readcsv("C1_M.csv"), readcsv("C2_M.csv"),
           readcsv("C3_M.csv"), readcsv("C4_M.csv"), readcsv("C5_M.csv")]

"""
Moulin discharge source for suite B.  Note that the moulin locations are in
the *.csv files in this folder.
"""
source_B_m(moulin_id, t, run) = moulins[run][moulin_id+1,end]
"Distributed source for suite B"
source_B(x,y,t,run) = source_A(x,y,t,1)

# C ################
const day = 24*60*60
const ra = [1/4, 1/2, 1, 2]
"""
Moulin discharge source for suite C.
"""
source_C_m(moulin_id,t,run) = max(0, source_B_m(moulin_id,t,5) * (1 - ra[run]*sin(2*pi*t/day)))

"""
Distributed source for suite C.
"""
source_C(x,y,t,run) = source_A(x,y,t,1)

# D ################
const year=31536000
const DT = [-4,-2,0,2,4] # delta T
const DDF = 0.01/day     # degree day factor m/k/s
const lr = -0.0075       # lapse rate in K/m

include("../topography/topos.jl") # need the surface elevation function

temp(t,run) = -16*cos(2*pi/year*t)- 5 + DT[run]

"""
Distributed source for suite D
"""
function source_D(x,y,t,run)
    z_s = sqrt_surface(x,y)
    return max(0, (z_s*lr+temp(t,run))*DDF) + source_A(x,y,t,1)
end
source_D_m(moulin_id,t,run) = 0

# E ################
"""
Steady distributed source of Suite E, in m/s.
"""
source_E(x,y,t,run) = 2 * source_A(x,y,t,6)
source_E_m(moulin_id,t,run) = 0

# F ################
const DT_valley = [-6,-3,0,3,6] # delta T
temp_valley(t,run) = -16*cos(2*pi/year*t)- 5 + DT_valley[run]

"""
Distributed source for suite F
"""
function source_F(x,y,t,run)
    z_s = valley_surface(x,y)
    return max(0, (z_s*lr+temp(t,run))*DDF) + source_A(x,y,t,1)
end
source_F_m(moulin_id,t,run) = 0

##########
# Plotting
##########
# set to `true` to enable plotting.  Requires the package PyPlot,
# install with `Pkg.add("PyPlot").
function plotsrc(source, source_m, x, y, ts, runs, tit)
    f = figure()
    subplot(2,1,1)
    title("Suite $tit")
    hold(true)
    for r in runs
        plot(ts/day, [source(x,y,t,r) for t in ts]*day, label="$r")
    end
    xlabel("Time (days)")
    ylabel("Dist source (m/day)")
    subplot(2,1,2)
    hold(true)
    moulin_id = 0
    for r in runs
        plot(ts/day, [source_m(moulin_id,t,r) for t in ts], label="$r")
    end
    xlabel("Time (days)")
    ylabel("Moulin source (m^3/s)")
    f
end
if true
    eval(:(using PyPlot))
    tday = linspace(0,day,500)
    tyear = linspace(0,year,500)
    plotsrc(source_A,source_A_m, 0,0, tday, 1:5, "A")
    plotsrc(source_B,source_B_m, 0,0, tday, 1:5, "B")
    ylim([0,95])
    plotsrc(source_C,source_C_m, 0,0, tday, 1:4, "C")
    plotsrc(source_D,source_D_m, 0,0, tyear, 1:5, "D")
    plotsrc(source_E,source_E_m, 0,0, tday, 1:5, "E")
    plotsrc(source_F,source_F_m, 0,0, tyear, 1:5, "F")
end
