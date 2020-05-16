module SailingNavigation

export Position, lat, lon, GridPoint, TimeSlice, TimedPath, closestpoint, surround
export RoutingProblem, minimumtimeroute

using GeometryTypes
using SailingPolars

# NB: This uses latitude (often considered to be y) first then longitude (often considered to be x).
# This latitude, then longitude ordering is as per ISO 6709 (en.wikipedia.org/wiki/ISO_6709)

# Position is a Float32 2-tuple of latitude and longitude in degrees
Position = Point2f0

# latitude from Position
lat(p::Position) = p[1]

# longitude from Position
lon(p::Position) = p[2]

# A GridPoint is a Position with the SurfaceParameters of wind and current at the Position
mutable struct GridPoint
    pt::Position
    sp::SurfaceParameters
end

"""
    TimeSlice

A TimeSlice is a matrix of GridPoints, each Position point with their SurfaceParameters
A Vector of TimeSlice can give the surface characteristics for an ocean region over time.
"""
TimeSlice = Matrix{GridPoint}

"""
    mutable struct RoutingProblem

timeinterval: the seconds duration for each TimeSlice
timeframe: a vector of sequential timeslices for the ocean region
obstacleindices: the Cartesian indices in each timeslice for
    obstacles, such as land or shoals, where the ship may not go
startindex: the timeslice position for time of starting
start: starting location on grid of GridPoints
finish: destination / finish location on grid of GridPoints
allowrepeatvisits: whether the vessel may overlap its prior path, usually false
"""
mutable struct RoutingProblem
    timeinterval::Float64 # minutes between timeframe slices
    timeframe::Vector{TimeSlice}
    obstacleindices::Vector{Point2{Int}}
    startindex::Int
    start::Point2{Int}
    finish::Point2{Int}
    allowrepeatvisits::Bool
end

"""
    mutable struct TimedPath

duration: minutes total to travel the path
path: vector of Cartesian indices of points in grid for path to travel
"""
mutable struct TimedPath
    duration::Float64
    path::Vector{Point2{Int}}
end

"""
    closestpoint(p, mat)

Get the closest GridPoint in matrix mat to a given position p.
p: Cartesian indices of a Position (latitude, longitude in degrees) in grid of GridPoints
mat: matrix of Gridpoints
"""
closestpoint(p, mat) = findmin(gp -> haversine(p[1], p[2], gp.pt[1], gp.pt[2])[1], mat)[2]

function surround(p, mat, excluded)
    neighbors = Point2{Int}[(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    (xmax, ymax) = size(mat)
    return filter(q -> 1 <= q[1] <= xmax && 1 <= q[2] <= ymax && !(q in excluded),
        [x + p for x in neighbors])
end

"""
    function minimumtimeroute(rp::RoutingProblem, sp::SailingPolar, verbose=false)

Get the route (as a TimedPath) that minimizes time from start to finish for a given
RoutingProblem (sea parameters) given sailing polar data (ship parameters).
"""
function minimumtimeroute(rp::RoutingProblem, sp::SailingPolar, verbose=false)
    timedpaths = [TimedPath(0.0, [rp.start])]
    completed, mintime, minpath = false, 1000.0, TimedPath(1000.0, [])
    for i in 1:1000
        newpaths = TimedPath[]
        verbose && println("Checking ", length(timedpaths), " paths of length ",
            length(timedpaths[1].path) - 1)
        for tpath in timedpaths
            if tpath.path[end] == rp.finish
                completed = true
                push!(newpaths, tpath)
            else
                p1 = tpath.path[end]
                slice = rp.timeframe[div(Int(round(tpath.duration)),
                                     Int(round(rp.timeinterval))) %
                                     length(rp.timeframe) + 1]
                for p2 in surround([p1[1], p1[2]], slice, rp.obstacleindices)
                    !rp.allowrepeatvisits && p2 in tpath.path && continue
                    gp1, gp2 = slice[p1[1], p1[2]], slice[p2[1], p2[2]]
                    lat1, lon1, lat2, lon2 = gp1.pt[1], gp1.pt[2], gp2.pt[1], gp2.pt[2]
                    t = sailsegmenttime(sp, gp1.sp, lat1, lon1, lat2, lon2)
                    path = deepcopy(tpath.path)
                    push!(path, p2)
                    push!(newpaths, TimedPath(tpath.duration + t, path))
                end
            end
        end
        timedpaths = unique(newpaths)
        if completed
            mindur = minimum(map(x -> x.duration, timedpaths))
            finished = filter(x -> x.path[end] == rp.finish, timedpaths)
            minfindur, idx = findmin(map(x -> x.duration, finished))
            verbose && println("Current finished minimum: ", finished[idx], ", others $mindur")
            if mindur == minfindur
                minpath = finished[idx]
                break
            end
        end
    end
    return minpath
end

end # module
