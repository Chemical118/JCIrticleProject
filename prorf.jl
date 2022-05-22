module prorf
    using ShapML, DataFrames, DecisionTree
    using Random, Statistics, Printf, PyPlot, StatsBase
    using FASTX, BioAlignments

    abstract type AbstractRF end
    abstract type AbstractRFI <: AbstractRF end

    struct RF <: AbstractRF
        fasta_loc::String
        data_loc::String
        amino_loc::Int
    end

    struct RFI <: AbstractRFI
        fasta_loc::String
        data_loc::String
        amino_loc::Int
        nfeat::StepRange{Int64, Int64}
        ntree::StepRange{Int64, Int64}
    end

    function RF(fasta_loc::String, data_loc::String)
        return RF(fasta_loc, data_loc, 0)
    end

    function RFI(fasta_loc::String, data_loc::String, nfeat::StepRange{Int64, Int64}, ntree::StepRange{Int64, Int64})
        return RFI(fasta_loc, data_loc, 0, nfeat, ntree)
    end

    function blosum_number(blosum_num::Int)
        if blosum_num == 45
            return BLOSUM45
        elseif blosum_num == 50
            return BLOSUM50
        elseif blosum_num == 62
            return BLOSUM62
        elseif blosum_num == 80
            return BLOSUM80
        elseif blosum_num == 90
            return BLOSUM90
        else
            error(@sprintf "There are no Matrix such as BLOSUM%d" blosum_num)
        end
    end

    function _find_key(d::Dict{Char, Int}, tar::Int)
        for k in keys(d)
            if d[k] == tar
                return k
            end
        end
    end

    function get_data(s::AbstractRF, ami_arr::Int, excel_loc::Char; norm::Bool=false, blosum::Int=62)
        data_len, loc_dict_vector, seq_matrix = _location_data(s.fasta_loc)
        blo = blosum_number(blosum)
        for (ind, (dict, col)) in enumerate(zip(loc_dict_vector, eachcol(seq_matrix)))
            max_num = maximum(values(dict))
            max_amino = _find_key(dict, max_num)
            if '-' ∉ keys(dict) && data_len - max_num ≥ ami_arr
                println(ind)
            end

        end
        # return Matrix{Float64}(A), Vector{Float64}(B), Vector{Tuple{Int, Char}}([(C[i, 1], only(C[i, 2])) for i = 1:size(C, 1)])
    end

    function get_amino_loc(s::AbstractRF, loc::Vector{Tuple{Int, Char}})
        return [string(i[1] + s.amino_loc) for i in loc]
    end
    
    function view_mutation(s::AbstractRF)
        _view_mutation(s.fasta_loc)
    end

    function view_mutation(fasta_loc::String)
        _view_mutation(fasta_loc)
    end

    function _view_mutation(fasta_loc::String)
        data_len, loc_dict_vector, _ = _location_data(fasta_loc)
        loc_vector = zeros(Int, data_len)
        loc_hist_vector = Vector()
        for dict in loc_dict_vector
            max_value = maximum(values(dict))
            if max_value ≠ data_len
                loc_vector[data_len - max_value] += 1
                push!(loc_hist_vector, data_len - max_value)
            end
        end
        loc_cum_vector = cumsum(loc_vector[end:-1:1])[end:-1:1]
        last_ind = findfirst(isequal(0), loc_cum_vector)
        loc_vector = loc_vector[1:last_ind]
        loc_cum_vector = loc_cum_vector[1:last_ind]
        plot(collect(1:last_ind), loc_cum_vector)
        xlabel("Number of mutaion")
        ylabel("Number of total amino location")
        PyPlot.title("Mutaion location stacked graph")
        xticks(collect(1:last_ind), [i % 5 == 1 ? i : "" for i in 1:last_ind])
        display(gcf())
        close("all")
        hist(loc_hist_vector, bins=20)
        xlim(1, last_ind - 1)
        xticks(collect(1:last_ind - 1), [i % 10 == 1 ? i : "" for i in 1:last_ind - 1])
        xlabel("Number of mutaion")
        ylabel("Number of amino location")
        PyPlot.title("Mutaion location histogram")
        display(gcf())
        close("all")
    end

    function _location_data(fasta_loc::String)
        reader = open(FASTA.Reader, fasta_loc)
        seq_vector = [collect(FASTA.sequence(String, record)) for record in reader]
        if length(Set(map(length, seq_vector))) ≠ 1
            error(@sprintf "%s is not aligned, Please align your data" fasta_loc)
        end
        seq_matrix = permutedims(hcat(seq_vector...))
        loc_vector = Vector{Dict{Char, Int}}()
        for col in eachcol(seq_matrix)
            loc_dict = Dict{Char, Int}()
            for val in col
                loc_dict[val] = get(loc_dict, val, 0) + 1
            end
            push!(loc_vector, loc_dict)
        end
        return size(seq_matrix)[1], loc_vector, seq_matrix
    end

    function get_reg_importance(s::AbstractRF, x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int; 
        val_mode::Bool=false, split_size::Float64=0.3, show_number::Int=20, imp_iter::Int=60,
        data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)), imp_state::Int64=rand(0:typemax(Int64)))
        
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        regr = RandomForestRegressor(n_trees=tree, n_subfeatures=feet, min_samples_leaf=1, rng=learn_state)
        DecisionTree.fit!(regr, x_train, y_train)

        if val_mode == false
            predict_test = DecisionTree.predict(regr, x_test)
            predict_train = DecisionTree.predict(regr, x_train)
            scatter(vcat(y_train, y_test), vcat(predict_train, predict_test), color=vcat(["orange" for i = 1:length(y_train)], ["blue" for i = 1:length(y_test)]))
            PyPlot.title("Random Forest Regression Result")
            xlabel("True Values")
            ylabel("Predictions")
            axis("equal")
            axis("square")
            xlim(0, xlim()[2])
            ylim(0, ylim()[2])
            plot([-1000, 1000], [-1000, 1000], color="black")
            display(gcf())
            close("all")
            @printf "NRMSE : %.6f\n" nrmse(predict_test, y_test)
        end
        return regr, _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), imp_iter, seed=imp_state, show_number=show_number, val_mode=val_mode)
    end

    function iter_get_reg_importance(s::AbstractRF, x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int, iter::Int;
        val_mode::Bool=false, split_size::Float64=0.3, show_number::Int=20, imp_iter::Int=60,
        data_state::Int64=rand(0:typemax(Int64)), imp_state::Int64=rand(0:typemax(Int64)))
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        f = zeros(length(loc), iter)
        n = zeros(iter)
        loc_list = get_amino_loc(s, loc)
        for i in 1:iter
            f[:, i], n[i] = _iter_get_reg_importance(x, x_train, x_test, y_train, y_test, loc_list, feet, tree, imp_iter, imp_state)
        end
        mf, sf = mean(f, dims=2)[:, 1], std(f, dims=2)[:, 1]

        if val_mode == false
            _iter_view_importance(mf, sf, loc_list, show_number=show_number)
            @printf "NRMSE : %.6f\n" mean(n)
        end

        return mf, sf
    end

    function _iter_get_reg_importance(x::Matrix{Float64}, x_train::Matrix{Float64}, x_test::Matrix{Float64}, y_train::Vector{Float64}, y_test::Vector{Float64}, loc::Vector{String}, feet::Int, tree::Int, imp_iter::Int, imp_state::Int64)
        regr = RandomForestRegressor(n_trees=tree, n_subfeatures=feet, min_samples_leaf=1)
        DecisionTree.fit!(regr, x_train, y_train)
        return _rf_importance(regr, DataFrame(x, loc), imp_iter, seed=imp_state, val_mode=true) , nrmse(regr, x_test, y_test)
    end
    
    function _iter_view_importance(fe::Vector{Float64}, err::Vector{Float64}, loc::Vector{String}; show_number::Int=20)
        norm_val = maximum(fe)
        fe /= norm_val
        err /= norm_val
        sorted_idx = sortperm(fe, rev=true)
        bar_pos = [length(sorted_idx):-1:1;] .- 0.5
        barh(bar_pos[1:show_number], fe[sorted_idx][1:show_number], xerr=err[sorted_idx][1:show_number], align="center", capsize=2)
        yticks(bar_pos[1:show_number], loc[sorted_idx][1:show_number])
        xlabel("Feature Importance")
        ylabel("Amino acid Location")
        PyPlot.title("Relative Mean Absolute Shapley Value")
        display(gcf())
        close("all")
    end

    function iter_view_importance(s::AbstractRF, loc::Vector{Tuple{Int, Char}}, fe::Vector{Float64}, err::Vector{Float64}; show_number::Int=20)
        _iter_view_importance(fe, err, get_amino_loc(s, loc), show_number=show_number)
    end

    function get_reg_value(si::AbstractRFI, x::Matrix{Float64}, y::Vector{Float64};
        val_mode::Bool=false, split_size::Float64=0.3, data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)))
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        z = zeros(Float64, length(si.nfeat), length(si.ntree))
        task = [(i[1], j[1], i[2], j[2]) for i in enumerate(si.nfeat), j in enumerate(si.ntree)]
        
        Threads.@threads for (i, j, feat, tree) in task
            z[i,  j] = _get_reg_value(x_train, x_test, y_train, y_test, feat, tree, learn_state)
        end

        if val_mode == false
            view_reg3d(si, z, title="NRMSE value")
        end
        return z
    end

    function get_reg_value_loc(si::AbstractRFI, z::Matrix{Float64})
        i, j = Tuple(findmin(z)[2])
        return collect(si.nfeat)[i], collect(si.ntree)[j] 
    end

    function _get_reg_value(x_train::Matrix{Float64}, x_test::Matrix{Float64}, y_train::Vector{Float64}, y_test::Vector{Float64}, feat::Int, tree::Int, learn_state::Int64)
        regr = RandomForestRegressor(n_trees=tree, n_subfeatures=feat, min_samples_leaf=1, rng=learn_state)
        DecisionTree.fit!(regr, x_train, y_train)
        return nrmse(regr, x_test, y_test)
    end

    function view_reg3d(s::AbstractRFI, z::Matrix{Float64}; title=nothing, elev=nothing, azim=nothing, scale::Int=2)
        nfeat_list = [s.nfeat;]' .* ones(length(s.ntree))
        ntree_list = ones(length(s.nfeat))' .* [s.ntree;]
        fig = figure()
        ax = fig.add_subplot(projection="3d")
        ax.view_init(elev=elev, azim=azim)
        xlabel("Numer of subfeatures")
        ylabel("Number of trees")
        mp = 10^scale
        ax.set_zlim(floor(minimum(z)*mp)/mp, ceil(maximum(z)*mp)/mp)
        ax.zaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
        surf = plot_surface(nfeat_list, ntree_list, z', cmap=ColorMap("coolwarm"), linewidth=0)
        colorbar(surf, shrink=0.7, aspect=15, pad=0.1, ax=ax)
        PyPlot.title(title)
        display(gcf())
        close("all")
    end

    function iter_get_reg_value(si::AbstractRFI, x::Matrix{Float64}, y::Vector{Float64}, iter::Int;
        val_mode::Bool=false, split_size::Float64=0.3, learn_state::Int64=rand(0:typemax(Int64)))
        z = zeros(length(si.nfeat), length(si.ntree), iter)
        for i = 1:iter
            z[:, :, i] = get_reg_value(si, x, y, val_mode=true, split_size=split_size, learn_state=learn_state)
        end

        vz, sz = mean(z, dims=3)[:, :, 1], std(z, dims=3)[:, :, 1]

        if val_mode == false
            view_reg3d(si, vz, title="NRMSE value", scale=2)
            view_reg3d(si, sz, title="NRMSE SD value", scale=3)
        end

        return vz, sz
    end

    function rf_importance(s::AbstractRF, regr::RandomForestRegressor, x::Matrix{Float64}, loc::Vector{Tuple{Int, Char}}, iter::Int=60; seed::Int64=rand(0:typemax(Int64)), val_mode::Bool=false)
        return _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), iter, seed=seed, val_mode=val_mode)
    end

    function _rf_dfpredict(regr::RandomForestRegressor, x::DataFrame)
        return DataFrame(y_pred = DecisionTree.predict(regr, Matrix{Float64}(x)))
    end

    function _rf_importance(regr::RandomForestRegressor, dx::DataFrame, iter::Int=60; 
                           seed::Int64=rand(0:typemax(Int64)), show_number::Int=20, val_mode::Bool=false)
        data_shap = ShapML.shap(explain = dx,
                        model = regr,
                        predict_function = _rf_dfpredict,
                        sample_size = iter,
                        seed = seed)
        data_plot = combine(groupby(data_shap, :feature_name), :shap_effect => x -> mean(abs.(x)))
        baseline = data_shap.intercept[1]
        feature_importance = data_plot[!, :shap_effect_function] / baseline
        if val_mode == false
            _view_importance(feature_importance, data_plot[!, :feature_name], baseline, show_number=show_number)
        end
        return feature_importance
    end

    function _view_importance(fe::Vector{Float64}, get_loc::Vector{String}; show_number::Int=20)
        fe /= maximum(fe)
        sorted_idx = sortperm(fe, rev=true)
        bar_pos = [length(sorted_idx):-1:1;] .- 0.5
        barh(bar_pos[1:show_number], fe[sorted_idx][1:show_number], align="center")
        yticks(bar_pos[1:show_number], get_loc[sorted_idx][1:show_number])
        xlabel("Feature Importance")
        ylabel("Amino acid Location")
        PyPlot.title("Relative Mean Absolute Shapley Value")
        display(gcf())
        close("all")
    end

    function _view_importance(fe::Vector{Float64}, get_loc::Vector{String}, baseline::Float64; show_number::Int=20)
        sorted_idx = sortperm(fe, rev=true)
        bar_pos = [length(sorted_idx):-1:1;] .- 0.5
        barh(bar_pos[1:show_number], fe[sorted_idx][1:show_number], align="center")
        yticks(bar_pos[1:show_number], get_loc[sorted_idx][1:show_number])
        xlabel(@sprintf "|Shapley effect| (baseline = %.2f)" baseline)
        ylabel("Amino acid Location")
        PyPlot.title("Feature Importance - Mean Absolute Shapley Value")
        display(gcf())
        close("all")
    end

    function view_importance(s::AbstractRF, loc::Vector{Tuple{Int, Char}}, fe::Vector{Float64}; show_number::Int=20)
        _view_importance(fe, get_amino_loc(s, loc), show_number=show_number)
    end

    function train_test_split(x::Matrix{Float64}, y::Vector{Float64}; test_size::Float64=0.3, random_state::Int64=rand(0:typemax(Int64)))
        # https://discourse.julialang.org/t/simple-tool-for-train-test-split/473/3
        n = length(y)
        idx = shuffle(MersenneTwister(random_state), 1:n)
        train_idx = view(idx, (floor(Int, test_size*n)+1):n)
        test_idx = view(idx, 1:floor(Int, test_size*n))
        x[train_idx,:], x[test_idx,:], y[train_idx], y[test_idx]
    end

    function nrmse(pre::Vector{Float64}, tru::Vector{Float64})
        return L2dist(pre, tru) / (maximum(tru) - minimum(tru)) / length(tru)^0.5
    end

    function nrmse(regr::RandomForestRegressor, test::Matrix{Float64}, tru::Vector{Float64})
        return L2dist(DecisionTree.predict(regr, test), tru) / (maximum(tru) - minimum(tru)) / length(tru)^0.5
    end

    function view_sequence(s::AbstractRF, typ::String="fasta", fontsize::Int=9, plot_width::Int=800)
        RFF.view_sequence(s.fasta_loc, s.amino_loc, typ=typ, fontsize=string(fontsize) * "pt", plot_width=plot_width)
    end
end