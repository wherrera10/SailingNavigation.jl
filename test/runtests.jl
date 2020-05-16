using GeometryTypes
using SailingPolars
using SailingNavigation
using Test

#=
The data is selected so the best time path is slightly longer than the
shortest length path. The forbidden regions are x, representing land or reef.
The allowed sailing points are . and start and finish are S and F.

x  .  .  F  .  .  x  .  x
.  .  .  .  .  .  .  x  x
x  .  .  x  x  x  .  .  .
.  .  x  x  x  x  .  x  x
x  .  .  .  x  x  .  x  .
x  .  .  .  x  x  .  x  .
.  .  .  .  x  .  .  x  .
x  .  .  .  .  .  .  x  .
.  .  .  S  .  .  .  .  .
=#
const forbidden = Point2{Int}.([
    [1, 8], [2, 1], [2, 8], [3, 5], [3, 8], [4, 1], [4, 5], [4, 6], [4, 8], [5, 1],
    [5, 5], [5, 6], [5, 8], [6, 3], [6, 4], [6, 5], [6, 6], [6, 8], [6, 9], [7, 1],
    [7, 4], [7, 5], [7, 6], [8, 8], [8, 9], [9, 1], [9, 7], [9, 9],
])

# Create regional wind patterns on the map.
function surfacebylongitude(lon)
    return lon < -155.03 ? SurfaceParameters(-5.0, 8, 150, 0.5) :
           lon < -155.99 ? SurfaceParameters(-90.0, 20, 150, 0.4) :
                           SurfaceParameters(180.0, 25, 150, 0.3)
end

# Vary wind speeds over time.
function mutatetimeslices!(slices)
    for (i, slice) in enumerate(slices), x in slice
        x.sp = SurfaceParameters(x.sp.winddeg, x.sp.windkts * (1 + 0.002 * i),
            x.sp.currentdeg, x.sp.currentkts)
    end
end


const startpos = Point2{Int}(1, 4)
const endpos = Point2{Int}(9, 4)
const pmat  = [Position(19.78 - 1/60 + i/60, -155.0 - 5/60 + j/60) for i in 0:8, j in 0:8]
const gpoints = map(pt -> GridPoint(pt, surfacebylongitude(lon(pt))), pmat)
const slices = [deepcopy(gpoints) for _ in 1:200]
mutatetimeslices!(slices)

const routeprob = RoutingProblem(10.0, slices, forbidden, 1, startpos, endpos, false)
const filename = "polar.csv"
const sp = getpolardata(filename)
const tp = minimumtimeroute(routeprob, sp)

@test tp.duration ≈ 4.0 * 60 + 43.69787966871
@test tp.path == Point{2,Int64}[[1, 4], [1, 5], [2, 6], [3, 7], [4, 7], [5, 7], [6, 7], [7, 7], [8, 6], [8, 5], [9, 4]]
