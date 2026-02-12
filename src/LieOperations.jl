export 
  generate_lyndon_words, 
  dim_free_nil_lie_alg 

function nextword(k::Int, w::Vector{Int}, alphabet::Vector{Int})
    # repeat w enough times and truncate to length n
    reps = (k ÷ length(w)) + 1
    x = repeat(w, reps)[1:k]
    # remove trailing maximal letters
    while !isempty(x) && x[end] == alphabet[end]
        pop!(x)
    end
    if !isempty(x)
        last_char = x[end]
        next_char_index = findfirst(==(last_char), alphabet) + 1
        x[end] = alphabet[next_char_index]
    end
    return x
end

function generate_lyndon_words(k::Int, alphabet::Vector{Int})
    lwords = Vector{Vector{Int}}()
    w = [alphabet[1]]
    while length(w) <= k
        push!(lwords, copy(w))
        w = nextword(k, w, alphabet)
        isempty(w) && break
    end
    return lwords
end

function dim_free_nil_lie_alg(d::Int,k::Int)
  return sum([div(moebius_mu(a)*d^(div(ell,a)),ell) for ell in (1:k) for a in (1:k) if divides(ell,a)[1]])
end
