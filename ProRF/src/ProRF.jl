module ProRF

using ShapML, DataFrames, DecisionTree
using PyCall, Random, Statistics, Printf, PyPlot, StatsBase
using FASTX, BioAlignments, XLSX, Phylo, AxisArrays

# export RF, RFI
# export blosum, get_data, get_amino_loc, view_mutation, view_reg3d, view_importance, view_sequence
# export train_test_split, nrmse, install_python_dependency
# export get_reg_importance, iter_get_reg_importance, iter_view_importance
# export get_reg_value, get_reg_value_loc, iter_get_reg_value, rf_importance
# export data_preprocess_fill, data_preprocess_index

# Struct defination

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
    return RF(fasta_loc, data_loc, 1)
end

function RFI(fasta_loc::String, data_loc::String, nfeat::StepRange{Int64, Int64}, ntree::StepRange{Int64, Int64})
    return RFI(fasta_loc, data_loc, 1, nfeat, ntree)
end

# Normal function

function blosum_matrix(blosum_num::Int)
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

function nrmse(pre::Vector{Float64}, tru::Vector{Float64})
    return L2dist(pre, tru) / (maximum(tru) - minimum(tru)) / length(tru)^0.5
end

function nrmse(regr::RandomForestRegressor, test::Matrix{Float64}, tru::Vector{Float64})
    return L2dist(DecisionTree.predict(regr, test), tru) / (maximum(tru) - minimum(tru)) / length(tru)^0.5
end

function predict_data(regr::RandomForestRegressor, test::Matrix{Float64})
    return DecisionTree.predict(regr, test)
end

function predict_data(regr::RandomForestRegressor, loc::Vector{Tuple{Int, Char}}, seq_vector::Vector{String}; blosum::Int=62)
    test_vector = Vector{Vector{Float64}}()
    blo = blosum_matrix(blosum)
    for seq in seq_vector
        push!(test_vector, [blo[tar, s] for ((_, tar), s) in zip(loc, seq)])
    end
    return DecisionTree.predict(regr, hcat(test_vector...)')
end

function view_sequence(s::AbstractRF; fontsize::Int=9, plot_width::Int=800, val_mode=false)
    seq_vector = [FASTA.sequence(String, record) for record in open(FASTA.Reader, s.fasta_loc)]
    id_vector = [FASTA.identifier(record) for record in open(FASTA.Reader, s.fasta_loc)]
    _view_sequence(s.fasta_loc, seq_vector, id_vector, s.amino_loc, fontsize * "pt", plot_width, val_mode=val_mode)
end

function view_sequence(fasta_loc::String, amino_loc::Int=0; fontsize::Int=9, plot_width::Int=800, val_mode=false)
    seq_vector = [FASTA.sequence(String, record) for record in open(FASTA.Reader, fasta_loc)]
    id_vector = [FASTA.identifier(record) for record in open(FASTA.Reader, fasta_loc)]
    _view_sequence(fasta_loc, seq_vector, id_vector, amino_loc, fontsize, plot_width, val_mode=val_mode)
end

function _view_sequence(fasta_loc::String, seq_vector::Vector{String}, id_vector::Vector{String}, amino_loc::Int=0, fontsize::Int=9, plot_width::Int=800; save_view::Bool=true)
    py"_view_sequence"(fasta_loc, seq_vector, id_vector, amino_loc, fontsize=string(fontsize) * "pt", plot_width=plot_width, save_view=save_view)
end

function train_test_split(x::Matrix{Float64}, y::Vector{Float64}; test_size::Float64=0.3, random_state::Int64=rand(0:typemax(Int64)))
    # https://discourse.julialang.org/t/simple-tool-for-train-test-split/473/3
    n = length(y)
    idx = shuffle(MersenneTwister(random_state), 1:n)
    train_idx = view(idx, (floor(Int, test_size*n)+1):n)
    test_idx = view(idx, 1:floor(Int, test_size*n))
    x[train_idx,:], x[test_idx,:], y[train_idx], y[test_idx]
end

function install_python_dependency()
    py"_python_install"("numpy")
    py"_python_install"("matplotlib")
    py"_python_install"("bokeh")
end

function _location_data(fasta_loc::String)
    seq_vector = [collect(FASTA.sequence(String, record)) for record in open(FASTA.Reader, fasta_loc)]
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

function _find_key(d::Dict{Char, Int}, tar::Int)
    for k in keys(d)
        if d[k] == tar
            return k
        end
    end
end

function _min_max_norm(ve::Vector{Float64})
    mi = minimum(ve)
    ma = maximum(ve)
    return [(i - mi) / (ma - mi) for i in ve]
end

function __init__()
    py"""
    def _python_install(package):
        import subprocess, sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

    def _view_sequence(floc, seqs, ids, loc, fontsize="9pt", plot_width=800, val_mode=False, save_view=True):
        # Bokeh sequence alignment view
        # https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner
        # Edit by Chemical118
        from bokeh.plotting import figure, output_file, show
        from bokeh.models import ColumnDataSource, Range1d
        from bokeh.models.glyphs import Text, Rect
        from bokeh.layouts import gridplot
        from bokeh.core.properties import value

        import numpy as np

        clrs = {'E': 'red', 'D': 'red', 'P': 'orange', 'A': 'orange', 'V': 'orange', 'H': 'orange', 'M': 'orange',
                'L': 'orange', 'I': 'orange', 'G': 'orange', 'K': 'blue', 'R': 'blue', 'N': 'green', 'C': 'green',
                'T': 'green', 'Q': 'green', 'S': 'green', 'F': 'yellow', 'Y': 'yellow', 'W': 'yellow', '-': 'white',
                'X': 'white', '.': 'white'}
        
        # make sequence and id lists from the aln object
        text = [it for s in list(seqs) for it in s]
        colors = [clrs[it] for it in text]
        n = len(seqs[0])
        s = len(seqs)
        # var = .4
        x = np.arange(1 + loc, n + 1 + loc)
        y = np.arange(0, s, 1)
        # creates a 2D grid of coords from the 1D arrays
        xx, yy = np.meshgrid(x, y)
        # flattens the arrays
        gx = xx.ravel()
        gy = yy.flatten()
        # use recty for rect coords with an offset
        recty = gy + .5
        # var = 1 / s
        # now we can create the ColumnDataSource with all the arrays
        source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
        plot_height = len(seqs) * 15 + 50
        x_range = Range1d(loc, n + loc + 1, bounds='auto')
        if n > 99:
            viewlen = 100
        else:
            viewlen = n + 1
        # view_range is for the close up view
        view_range = (0 + loc, viewlen + loc)
        tools = "xpan, xwheel_zoom, reset, save"

        # entire sequence view (no text, with zoom)
        p = figure(title=None, plot_width=plot_width, plot_height=50,
                x_range=x_range, y_range=(0, s), tools=tools,
                min_border=0, toolbar_location='below')
        rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                    line_color=None, fill_alpha=0.6)
        p.add_glyph(source, rects)
        p.yaxis.visible = False
        p.grid.visible = False

        # sequence text view with ability to scroll along x-axis
        p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
                    x_range=view_range, y_range=ids, tools="xpan,reset",
                    min_border=0, toolbar_location='below')  # , lod_factor=1)
        glyph = Text(x="x", y="y", text="text", text_align='center', text_color="black",
                    text_font=value("arial"), text_font_size=fontsize)
        rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                    line_color=None, fill_alpha=0.4)
        p1.add_glyph(source, glyph)
        p1.add_glyph(source, rects)

        p1.grid.visible = False
        p1.xaxis.major_label_text_font_style = "bold"
        p1.yaxis.minor_tick_line_width = 0
        p1.yaxis.major_tick_line_width = 0

        p = gridplot([[p], [p1]], toolbar_location='below')

        if save_view:
            output_file('Data/View/' + floc.split('/')[-1].split('.')[0] + '.html')
        
        show(p)
    """
end

# Data preprocess

function data_preprocess_index(in_fasta_loc::String; target_rate::Float64=0.3, val_mode::Bool=false)
    data_len, loc_dict_vector, _ = _location_data(in_fasta_loc)

    seq_vector = Vector{String}()
    id_vector = Vector{String}()

    for record in open(FASTA.Reader, in_fasta_loc)
        push!(id_vector, FASTA.identifier(record))
        push!(seq_vector, FASTA.sequence(String, record))
    end

    gap_fre_vector = [get(dict, '-', 0) / data_len for dict in loc_dict_vector]

    front_ind = findfirst(x -> x ≤ target_rate, gap_fre_vector)
    last_ind = findlast(x -> x ≤ target_rate, gap_fre_vector)
    seq_len = last_ind - front_ind + 1
    
    seq_vector = [seq[front_ind:last_ind] for seq in seq_vector]

    if val_mode == false
        if seq_len ≥ 100
            _view_sequence(in_fasta_loc, [seq[1:45] * '.' ^ 9 * seq[end-44:end] for seq in seq_vector], id_vector, save_view=false)
        else
            _view_sequence(in_fasta_loc, seq_vector, id_vector, save_view=false)
        end
    end

    return front_ind, last_ind
end

function data_preprocess_fill(front_ind::Int, last_ind::Int, in_fasta_loc::String, newick_loc::String, out_fasta_loc::String; val_mode::Bool=false)
    dis_matrix = distances(open(parsenewick, newick_loc))
    seq_vector = Vector{String}()
    id_vector = Vector{String}()
    id_dict = Dict{String, Int}()

    for (ind, record) in enumerate(open(FASTA.Reader, in_fasta_loc))
        record_id = FASTA.identifier(record)
        id_dict[record_id] = ind
        push!(id_vector, record_id)
        push!(seq_vector, FASTA.sequence(String, record)[front_ind:last_ind])
    end

    if length(Set(map(length, seq_vector))) ≠ 1
        error(@sprintf "%s is not aligned, Please align your data" fasta_loc)
    end

    seq_len = length(seq_vector[1])
    if Set(id_vector) ≠ Set(AxisArrays.axes(dis_matrix, 1))
        error("Make sure the data same name of tree and id of fasta file")
    end

    gap_ind_vector = Vector{Tuple{Int, Int}}()
    isgap_vector = Vector{Bool}()
    for seq in seq_vector
        front_gap_ind = findfirst(!isequal('-'), seq)
        last_gap_ind = findlast(!isequal('-'), seq)
        push!(isgap_vector, front_gap_ind > 2 || last_gap_ind < seq_len - 1)
        push!(gap_ind_vector, (front_gap_ind, last_gap_ind))
    end

    if any(isgap_vector)
        println("<Closest Target --> Main>")
    else
        println("There are no gap in data!")
        return
    end

    dis_matrix_id_vector = AxisArrays.axes(dis_matrix, Axis{1})
    dis_matrix_isgap_vector = [isgap_vector[id_dict[id]] for id in dis_matrix_id_vector]

    edit_seq_vector = Vector{String}()
    for (ind, (front_gap_ind, last_gap_ind)) in enumerate(gap_ind_vector)
        main_seq = seq_vector[ind]
        if isgap_vector[ind] == true
            nogap_dis_vector = map(i -> Float64(dis_matrix_isgap_vector[i[1]] == true ? 1 : i[2]), enumerate(dis_matrix[:, id_vector[ind]]))
            min_target_ind = id_dict[dis_matrix_id_vector[argmin(nogap_dis_vector)]]
            min_target_seq = seq_vector[min_target_ind]

            if front_gap_ind > 2
                main_seq = min_target_seq[1:front_gap_ind - 1] * main_seq[front_gap_ind:end]
            end

            if last_gap_ind < seq_len - 1
                main_seq = main_seq[1:last_gap_ind] * min_target_seq[last_gap_ind + 1:end]
            end
            @printf "%s --> %s\n" id_vector[min_target_ind] id_vector[ind]
        end
        push!(edit_seq_vector, main_seq)
    end

    open(FASTA.Writer, out_fasta_loc) do io
        for (seq, id) in zip(edit_seq_vector, id_vector)
            write(io, FASTA.Record(id, seq))
        end
    end

    if val_mode == false
        _view_sequence(out_fasta_loc, seq_vector, id_vector, save_view=false)
    end
end

# RF / RFI function

function get_data(s::AbstractRF, ami_arr::Int, excel_col::Char; norm::Bool=false, blosum::Int=62, sheet="Sheet1", title::Bool=true)
    _get_data(s, ami_arr, excel_col, norm, blosum, sheet, title)
end

function get_data(s::AbstractRF, excel_col::Char; norm::Bool=false, blosum::Int=62, sheet="Sheet1", title::Bool=true)
    _get_data(s, 1, excel_col, norm, blosum, sheet, title)
end

function _get_data(s::AbstractRF, ami_arr::Int, excel_col::Char, norm::Bool, blonum::Int, sheet::String, title::Bool)
    data_len, loc_dict_vector, seq_matrix = _location_data(s.fasta_loc)
    blo = blosum_matrix(blonum)
    x_col_vector = Vector{Vector{Float64}}()
    loc_vector = Vector{Tuple{Int, Char}}()
    for (ind, (dict, col)) in enumerate(zip(loc_dict_vector, eachcol(seq_matrix)))
        max_val = maximum(values(dict))
        max_amino = _find_key(dict, max_val)
        if '-' ∉ keys(dict) && ami_arr ≤ data_len - max_val 
            push!(x_col_vector, [blo[max_amino, i] for i in col])
            push!(loc_vector, (ind, max_amino))
            print(seq_matrix[1, ind])
        end
    end

    excel_data = DataFrame(XLSX.readtable(s.data_loc, sheet, infer_eltypes=title)...)
    excel_select_vector = excel_data[!, Int(excel_col) - Int('A') + 1]
    if norm
        excel_select_vector = _min_max_norm(excel_select_vector)
    end
    
    println(size(hcat(x_col_vector...)))
    x = Matrix{Float64}(hcat(x_col_vector...))
    y = Vector{Float64}(excel_select_vector)
    l = Vector{Tuple{Int, Char}}(loc_vector)
    return x, y, l
end

function get_amino_loc(s::AbstractRF, loc::Vector{Tuple{Int, Char}})
    return [string(i[1] + s.amino_loc - 1) for i in loc]
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

function rf_importance(s::AbstractRF, regr::RandomForestRegressor, x::Matrix{Float64}, loc::Vector{Tuple{Int, Char}}, iter::Int=60; seed::Int64=rand(0:typemax(Int64)), val_mode::Bool=false)
    return _rf_importance(regr, DataFrame(x, get_amino_loc(s, loc)), iter, seed=seed, val_mode=val_mode)
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

function _rf_dfpredict(regr::RandomForestRegressor, x::DataFrame)
    return DataFrame(y_pred = DecisionTree.predict(regr, Matrix{Float64}(x)))
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

function iter_view_importance(s::AbstractRF, loc::Vector{Tuple{Int, Char}}, fe::Vector{Float64}, err::Vector{Float64}; show_number::Int=20)
    _iter_view_importance(fe, err, get_amino_loc(s, loc), show_number=show_number)
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

# RFI function

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

function _get_reg_value(x_train::Matrix{Float64}, x_test::Matrix{Float64}, y_train::Vector{Float64}, y_test::Vector{Float64}, feat::Int, tree::Int, learn_state::Int64)
    regr = RandomForestRegressor(n_trees=tree, n_subfeatures=feat, min_samples_leaf=1, rng=learn_state)
    DecisionTree.fit!(regr, x_train, y_train)
    return nrmse(regr, x_test, y_test)
end

function get_reg_value_loc(si::AbstractRFI, z::Matrix{Float64})
    i, j = Tuple(findmin(z)[2])
    return collect(si.nfeat)[i], collect(si.ntree)[j] 
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

end