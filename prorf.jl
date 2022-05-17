module prorf
    using ShapML, DataFrames, DecisionTree, MLJ
    using PyCall, Random, Statistics, Printf
    using PyPlot
    
    pushfirst!(pyimport("sys")."path", "")
    RFC = pyimport("prorf.rfclass")
    RFF = pyimport("prorf.rfunction")

    struct RF
        amino_loc::Int
        fasta_loc::String
        data_loc::String
    end

    struct RFI
        rf::RF
        nfeat::StepRange{Int64, Int64}
        ntree::StepRange{Int64, Int64}
    end

    function RF(fasta_loc::String, data_loc::String)
        return RF(0, fasta_loc, data_loc)
    end

    function RFI(amino_loc::Int, fasta_loc::String, data_loc::String, nfeat::StepRange{Int64, Int64}, ntree::StepRange{Int64, Int64})
        return RFI(RF(amino_loc, fasta_loc, data_loc), nfeat, ntree)
    end

    function RFI(fasta_loc::String, data_loc::String, nfeat::StepRange{Int64, Int64}, ntree::StepRange{Int64, Int64})
        return RFI(RF(0, fasta_loc, data_loc), nfeat, ntree)
    end

    function get_data(s::RF, ami_arr::Int, excel_loc::Char, norm::Bool=false) where T<:RF
        R = RFC.RF(s.fasta_loc, s.data_loc, norm)
        A, B, C = R.get_data(ami_arr, string(excel_loc), norm)
        return Matrix{Float64}(A), Vector{Float64}(B), Vector{Tuple{Int, Char}}([(C[i, 1], only(C[i, 2])) for i = 1:size(C, 1)])
    end

    function get_data(si::RFI, ami_arr::Int, excel_loc::Char, norm::Bool=false) where T<:RF
        s = si.rf
        R = RFC.RF(s.fasta_loc, s.data_loc, norm)
        A, B, C = R.get_data(ami_arr, string(excel_loc), norm)
        return Matrix{Float64}(A), Vector{Float64}(B), Vector{Tuple{Int, Char}}([(C[i, 1], only(C[i, 2])) for i = 1:size(C, 1)])
    end

    function get_amino_loc(s::RF, loc::Vector{Tuple{Int, Char}})
        return [string(i[1] + s.amino_loc) for i in loc]
    end

    function get_amino_loc(si::RFI, loc::Vector{Tuple{Int, Char}})
        return [string(i[1] + si.rf.amino_loc) for i in loc]
    end

    function get_reg_importance(s::RF, x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int; 
        val_mode::Bool=false, split_size::Float64=0.3, show_number::Int=20, imp_iter::Int=60,
        data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)), imp_state::Int64=rand(0:typemax(Int64)))
        
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        regr = RandomForestRegressor(n_trees=tree, max_depth=feet, min_samples_leaf=1, rng=learn_state)
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
            @printf "NRMSE : %.6f" nrmse(predict_test, y_test)
        end
        return regr, _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), imp_iter, seed=imp_state, show_number=show_number, val_mode=val_mode)
    end

    function get_reg_importance(s::RF, x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int; 
        val_mode::Bool=false, split_size::Float64=0.3, show_number::Int=20, imp_iter::Int=60,
        data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)), imp_state::Int64=rand(0:typemax(Int64)))
        
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        regr = RandomForestRegressor(n_trees=tree, max_depth=feet, min_samples_leaf=1, rng=learn_state)
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
            @printf "NRMSE : %.6f" nrmse(predict_test, y_test)
        end
        return regr, _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), imp_iter, seed=imp_state, show_number=show_number, val_mode=val_mode)
    end

    function get_reg_importance(si::RFI, x::Matrix{Float64}, y::Vector{Float64}, loc::Vector{Tuple{Int, Char}}, feet::Int, tree::Int; 
        val_mode::Bool=false, split_size::Float64=0.3, show_number::Int=20, imp_iter::Int=60,
        data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)), imp_state::Int64=rand(0:typemax(Int64)))
        s = si.rf
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        regr = RandomForestRegressor(n_trees=tree, max_depth=feet, min_samples_leaf=1, rng=learn_state)
        DecisionTree.fit!(regr, x_train, y_train)

        if val_mode == false
            predict_test = DecisionTree.predict(regr, x_test)
            predict_train = DecisionTree.predict(regr, x_train)
            scatter(vcat(y_train, y_test), vcat(predict_train, predict_test), color=vcat(["orange" for i = 1:length(y_train)], ["blue" for i = 1:length(y_test)]))
            # PyPlot.title("Random Forest Regression Result")
            PyPlot.title(nothing)
            xlabel("True Values")
            ylabel("Predictions")
            axis("equal")
            axis("square")
            xlim(0, xlim()[2])
            ylim(0, ylim()[2])
            plot([-1000, 1000], [-1000, 1000], color="black")
            display(gcf())
            close("all")
            @printf "NRMSE : %.6f" nrmse(predict_test, y_test)
        end
        return regr, _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), imp_iter, seed=imp_state, show_number=show_number, val_mode=val_mode)::Vector{Float64}
    end

    function get_reg_value(si::RFI, x::Matrix{Float64}, y::Vector{Float64};
        val_mode::Bool=false, split_size::Float64=0.3, data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)))
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=split_size, random_state=data_state)
        z = zeros(Float64, length(si.nfeat), length(si.ntree))
        task = [(i[1], j[1], i[2], j[2]) for i in enumerate(si.nfeat), j in enumerate(si.ntree)]
        
        Threads.@threads for (i, j, feat, tree) in task
            regr = RandomForestRegressor(n_trees=tree, max_depth=feat, min_samples_leaf=1, rng=learn_state)
            DecisionTree.fit!(regr, x_train, y_train)
            z[i,  j] = nrmse(regr, x_test, y_test)
            # @printf "%d %d %d \n" Threads.threadid() feat tree
        end

        if val_mode == false
            view_reg3d(si, z)
        end

        return z
    end

    function view_reg3d(s::RFI, z::Matrix{Float64}; title="NRMSE value", c_format="%.3f", elev=nothing, azim=nothing, scale::Int=2)
        nfeat_list = [s.nfeat;]' .* ones(length(s.ntree))
        ntree_list = ones(length(s.nfeat))' .* [s.ntree;]
        fig = figure()
        ax = fig.add_subplot(projection="3d")
        ax.view_init(elev=elev, azim=azim)
        xlabel("Max tree depth")
        ylabel("Number of trees")
        mp = 10^scale
        ax.set_zlim(floor(minimum(z)*mp)/mp, ceil(maximum(z)*mp)/mp)
        ax.zaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
        surf = plot_surface(nfeat_list, ntree_list, z', cmap=ColorMap("coolwarm"), linewidth=0)
        colorbar(surf, shrink=0.7, aspect=15, pad=0.1, format=c_format, ax=ax)
        PyPlot.title(title)
        display(gcf())
        close("all")
    end

    # function iter_get_reg_value(si::RFI, x::Matrix{Float64}, y::Vector{Float64}, iter::Int;
    #     val_mode::Bool=false, split_size::Float64=0.3, data_state::Int64=rand(0:typemax(Int64)), learn_state::Int64=rand(0:typemax(Int64)))
        
        
    # end
    function rf_importance(s::RF, regr::RandomForestRegressor, x::Matrix{Float64}, loc::Vector{Tuple{Int, Char}}, iter::Int=60; seed::Int64=rand(0:typemax(Int64)), val_mode::Bool=false)
        return _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), iter, seed=seed, val_mode=val_mode)
    end

    function rf_importance(s::RFI, regr::RandomForestRegressor, x::Matrix{Float64}, loc::Vector{Tuple{Int, Char}}, iter::Int=60; seed::Int64=rand(0:typemax(Int64)), val_mode::Bool=false)
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

    function view_importance(s::RF, fe::Vector{Float64}, loc::Vector{Tuple{Int, Char}}; show_number::Int=20)
        _view_importance(fe, get_amino_loc(s, loc), show_number=show_number)
    end

    function view_importance(si::RFI, fe::Vector{Float64}, loc::Vector{Tuple{Int, Char}}; show_number::Int=20)
        _view_importance(fe, get_amino_loc(si.rf, loc), show_number=show_number)
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
        return rms(pre, tru)^0.5 / (maximum(tru) - minimum(tru))
    end

    function nrmse(regr::RandomForestRegressor, test::Matrix{Float64}, tru::Vector{Float64})
        pre = DecisionTree.predict(regr, test)
        return rms(pre, tru) / (maximum(tru) - minimum(tru))
    end

end