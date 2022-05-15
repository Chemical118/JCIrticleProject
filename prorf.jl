module prorf
    using PyCall
    using Random
    pushfirst!(pyimport("sys")."path", "")
    RFC = pyimport("prorf.rfclass")

    struct RF
        amino_loc::Int
        fasta_loc::String
        data_loc::String
    end

    function RF(fasta_loc::String, data_loc::String)
        return RF(0, fasta_loc, data_loc)
    end

    function get_data(s::RF, ami_arr::Int, excel_loc::Char, norm::Bool=false)
        R = RFC.RF(s.fasta_loc, s.data_loc, norm)
        A, B, C = R.get_data(ami_arr, string(excel_loc), norm)
        return Matrix{Float64}(A), Vector{Float64}(B), Vector{Tuple{Int, Char}}([(C[i, 1], only(C[i, 2])) for i = 1:size(C, 1)])
    end

    function get_reg_importance(x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int; 
                                split_size::Float64=0.3, )
        
    end

    function train_test_split(x::Matrix{Float64}, y::Vector{Float64}; test_size::Float64, random_state::Int)
        # https://discourse.julialang.org/t/simple-tool-for-train-test-split/473/3
        n = length(y)
        idx = shuffle(MersenneTwister(random_state), 1:n)
        train_idx = view(idx, (floor(Int, test_size*n)+1):n)
        test_idx = view(idx, 1:floor(Int, test_size*n))
        x[train_idx,:], x[test_idx,:], y[train_idx,:], y[test_idx,:]
    end
end