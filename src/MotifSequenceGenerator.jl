"""
    MotifSequenceGenerator
This module generates random sequences of motifs, under the constrain that the
sequence has some total length **exactly equal** to `q`.

The motifs are contained in a vector, and they can be *anything*.
The notion of "temporal length" is defined
based on the following two functions, which must be provided for the motifs:

* `limits(motif)` : Some function that given the `motif` it returns the
  `(start, end)` of the the motif in the same units as `q`. Defaults to
  `(1, length(motif))`. Notice that this function establishes a measure of
  temporal length, which simply is `end - start`.
* `translate(motif, t)` : Some function that given the `motif` it returns a *new*
  motif which is temporally translated by `t` (either negative or positive), with
  respect to the same units as `q`. Defaults to `motif .+ n`.

All main functionality is given by the function [`random_sequence`](@ref).

**At the moment this module has not been tested with non-integer lengths**.

TODO: For non-integer `q` ask a `δq` so that lengths up to `q ± δq` are accepted as
ok!
"""
module MotifSequenceGenerator

using Base.Iterators
using Combinatorics, Random

export random_sequence, all_possible_sums

struct DeadEndMotifs <: Exception
  tries::Int
  summands::Int
  tailcut::Int
end
Base.showerror(io::IO, e::DeadEndMotifs) = print(io,
"DeadEndMotifs: Couldn't find a proper sequence with $(e.tries) random tries, "*
"each with summands up to $(e.summands) (total tailcuts: $(e.tailcut)).")



default_limits(motif) = (1, length(motif))
default_translate(motif, t) = motif .+ t



"""
    random_sequence(motifs::Vector{M}, q [, limits, translate]; kwargs...)
Create a random sequence of motifs of type `M`, under the constraint that the
sequence has "length" **exactly** `q`. Return the sequence itself as well as the
sequence of indices of `motifs` used to create it.

"length" here means an abstracted "temporal length" defined by the struct `M`,
based on the `limits` and `translate` functions.
It does **not** refer to the amount of elements, although coincidentally this is the
default behavior. Please see [`MotifSequenceGenerator`](@ref) for defining `limits`
and `translate`.

## Description & Keywords
The algorithm works as follows: First a random sequence of motifs is created,
so that it has length of `q ≤ s ≤ q - maximum(motiflengths)`. The possible tries
of random sequences is set by the `tries` keyword (default 5).

For each random try, it is first check whether the sequence is already correct.
If not, the last entry of the sequence is dropped. Then, since the sequence is now
already smaller than `q`, all possible sums of `summands` out of the motif pool
are checked. If some combination of `summands` sums to exactly the difference, they are
added to the sequence. For multiple satisfactory combinations, a random one is picked.

If the random combination of `summands` does not fit, one more entry is dropped
from the sequence and the process is repeated.

The keyword `tailcut` limits how many times will an element
be dropped (default is `2`).
The keyword `summands` denotes how many possible combinations of sums to check.
Default is `3` which means that we check from 2 up to `3` summands (*all*
possible combinations of summands are checked, which means that the computation
skyrockets for large `summands`!).

If after going though all these combinations of possible sequences we do not find
a proper one, an error is thrown.
"""
function random_sequence(motifs::Vector{M}, q::Int,
    limits = default_limits, translate = default_translate;
    tries = 5, summands = 3, tailcut = 2) where {M}

    idxs = 1:length(motifs)
    motifs0, motiflens = _motifs_at_origin(motifs, limits, translate)

    q < minimum(motiflens) && throw(ArgumentError(
    "Minimum length of motifs is less than `q`. Impossible to make a sequence."
    ))

    worked = false; count = 0; seq = Int[]
    while worked == false
        count > tries && throw(DeadEndMotifs(tries, summands, tailcut))

        seq, seq_length = _random_sequence_try(motiflens, q)
        worked = _complete_sequence!(seq, motiflens, q, summands, tailcut)
        count += 1
    end

    return _instantiate_sequence(motifs0, motiflens, seq, translate), seq
end

"""
    _motifs_at_origin(motifs, limits, translate) -> (motifs0, motiflens)
Bring all motifs to the origin and compute the motif lengths.
"""
function _motifs_at_origin(motifs::Vector{M}, limits, translate) where M
    motifs0 = similar(motifs)
    motiflens = zeros(Int, length(motifs))
    for i in 1:length(motifs)
        start, fine = limits(motifs[i])
        motifs0[i] = start == 0 ? motifs[i] : translate(motifs[i], -start)
        motiflens[i] = fine - start
    end
    return motifs0, motiflens
end

"""
    _random_sequence_try(motiflens, q) -> seq, seq_length
Return a random sequence of motif indices
so that the total sequence is *guaranteed* to have total length of
`q ≤ s ≤ q - maximum(motiflens)`.
"""
function _random_sequence_try(motiflens, q)
    seq = Int[]; seq_length = 0; idxs = 1:length(motiflens)
    while seq_length < q
        i = rand(idxs)
        push!(seq, i)
        seq_length += motiflens[i]
    end
    return seq, seq_length
end


function _complete_sequence!(seq, motiflens, q, summands, tailcut)

    remainder = q - sum(motiflens[k] for k in seq)
    if remainder == 0
        # Case 0: The sequence is already exactly equal to q
        return true
    elseif remainder < 0 && -remainder ∈ motiflens
        # Case 1: There is an extra difference, which is an
        # exact length of some motif.
        # We find the possible motifs, pick a random one, and pick
        # a random position in the sequence that it exists.
        # Delete that entry of the sequence.
        mi = rand(findall((in)(-remainder), motiflens))
        possible = findall((in)(mi), seq)
        if !isempty(possible)
            deleteat!(seq, rand(possible))
            return true
        end
    else
        # Case 2: Recursive deletion of last entry of the sequence, and trying to
        # see if it can be completed with some combination of existing motifs
        tcut = 0
        while tcut < tailcut
            tcut += 1
            pop!(seq)
            isempty(seq) && return false
            remainder = q - sum(motiflens[k] for k in seq)
            if remainder ∈ motiflens
                mi = rand(findall(in(remainder), motiflens))
                push!(seq, mi)
                return true
            end
            for n in 2:summands
                everything = all_possible_sums(motiflens, n)
                sums = [e[1] for e in everything]
                if remainder ∈ sums
                    cases = findall(in(remainder), sums)
                    if !isempty(cases)
                        idxs_of_vals = shuffle!(everything[rand(cases)][2])
                        push!(seq, idxs_of_vals...)
                        return true
                    end
                end
            end
        end
    end
    return false
end

# Function provided by Mark Birtwistle in stackoverflow
"""
    all_possible_sums(summands, n)
Compute all possible sums from combining `n` elements from `summands`,
with repetition and using only unique combinations.

Return a vector of tuples: the first
entry of each tuple is the sum, while the second is the indices of summands
used to compute the sum.
"""
function all_possible_sums(summands::Vector{T}, n) where {T}
    m = length(summands)
    r = Vector{Tuple{T, Vector{Int}}}(undef, binomial(m + n - 1, n))
    s = with_replacement_combinations(eachindex(summands), n)
    for (j, i) in enumerate(s)
        r[j] = (sum(summands[j] for j in i), i)
    end
    return r
end

function _instantiate_sequence(motifs0::Vector{M}, motiflens, seq, translate) where M
    ret = M[]
    prev = 0
    for s in seq
        push!(ret, translate(motifs0[s], prev))
        prev += motiflens[s]
    end
    return ret
end


end
