module prorf
    using ShapML, DataFrames, DecisionTree
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

    function get_reg_importance(s::RF, x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int; 
                                val_mode::Bool=false, split_size::Float64=0.3, show_number::Int=20, data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)))
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        regr = RandomForestRegressor(n_trees=tree, max_depth=feet, min_samples_leaf=1, rng=learn_state)
        fit!(regr, x_train, y_train)
        dx = DataFrame(x, [string(i[1]) for i in loc])
        rf_importance(regr, dx, show_number=show_number, val_mode=val_mode)
    end

    function rf_dfpredict(regr::RandomForestRegressor, x::DataFrame)
        return DataFrame(y_pred = predict(regr, Matrix{Float64}(x)))
    end

    function rf_importance(regr::RandomForestRegressor, dx::DataFrame; show_number::Int=20, val_mode::Bool=false)
        data_shap = ShapML.shap(explain = dx,
                        model = regr,
                        predict_function = rf_dfpredict,
                        sample_size = 60,
                        seed = 1)
    end

    function train_test_split(x::Matrix{Float64}, y::Vector{Float64}; test_size::Float64=0.3, random_state::Int64=rand(0:typemax(Int64)))
        # https://discourse.julialang.org/t/simple-tool-for-train-test-split/473/3
        n = length(y)
        idx = shuffle(MersenneTwister(random_state), 1:n)
        train_idx = view(idx, (floor(Int, test_size*n)+1):n)
        test_idx = view(idx, 1:floor(Int, test_size*n))
        x[train_idx,:], x[test_idx,:], y[train_idx], y[test_idx]
    end
end