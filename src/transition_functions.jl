"""
    lerp(x, x0, x1, y0, y1)
Interpolate between two points (x0, y0) and (x1, y1)
```jldoctest;setup = :(using HallThruster: lerp)
julia> lerp(0.5, 0.0, 1.0, 0.0, 2.0)
1.0
"""
@inline function lerp(x, x0, x1, y0, y1)
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end

abstract type TransitionFunction end

struct SmoothIf <: TransitionFunction
    transition_length::Float64
end

(s::SmoothIf)(x, cutoff, y1, y2) = ((y2 - y1)*tanh((x-cutoff)/s.transition_length) + y1 + y2) / 2

struct QuadraticTransition <: TransitionFunction
    transition_length::Float64
    offset::Float64
end

function (q::QuadraticTransition)(x, cutoff, y1, y2)
    x′ = x - q.offset
    if x′ < cutoff
        return y1
    else
        return y1 + (y2 - y1) * (x′ - cutoff)^2 / ((x - cutoff)^2 + q.transition_length^2)
    end
end

struct LinearTransition <: TransitionFunction
    transition_length::Float64
    offset::Float64
end
function(ℓ::LinearTransition)(x, cutoff, y1, y2)
    x′ = x - ℓ.offset
    L = ℓ.transition_length
    x1 = cutoff - L/2
    x2 = cutoff + L/2
    if x′ < x1
        return y1
    elseif x′ > x2
        return y2
    else
        return lerp(x, x1, x2, y1, y2)
    end
end

struct StepFunction <: TransitionFunction end

(::StepFunction)(x, cutoff, y1, y2) = x ≤ cutoff ? y1 : y2
