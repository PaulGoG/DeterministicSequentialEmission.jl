#s-wave neutron strength function from RIPL3 vs different parametrizations

using CSV, DataFrames, Plots, ElectronDisplay

rawdatafile_name = "resonances0.dat"
rawdatafile_header = ["Z", "Symbol", "A", "Io", "Bn", "D0", "dD", "S0", "dS", "Gg", "dG", "Com"]
rawdatafile_firstline = 5
rawdatafile_delim = ' '
rawdatafile = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, 
                        skipto = rawdatafile_firstline, header = rawdatafile_header)
rawdatafile.S0 .*= 1e-4

RIPL3_data = DataFrame(A = rawdatafile.A[rawdatafile.S0 .> 0], 
                    Z = rawdatafile.Z[rawdatafile.S0 .> 0], 
                    S_0 = rawdatafile.S0[rawdatafile.S0 .> 0])

function Strength_function_new(A)
    if A <= 23
        return 2.5e-5
    elseif A > 23 && A <= 55
        return 1.29695e-5 *A - 2.7328e-4
    elseif A > 55 && A <= 69
        return -2.1429e-5 *A + 1.6186e-3
    elseif A > 69 && A <= 79
        return 1.4e-4
    elseif A > 79 && A <= 90
        return -7.8182e-6 *A + 7.5764e-4
    elseif A > 90 && A <= 110
        return 5.4e-5
    elseif A > 110 && A <= 120
        return -4.4e-6 *A + 5.38e-4 
    elseif A > 120 && A <= 135
        return 6.6667e-6 *A - 7.9e-4
    elseif A > 135 && A <= 146
        return 2.4545e-5 *A - 3.2036e-3
    elseif A > 146 && A <= 152
        return 3.8e-4
    elseif A > 152 && A <= 158
        return -3.8333e-5 *A + 6.2067e-3
    elseif A > 158 && A <= 173
        return 1.5e-4
    elseif A > 173 && A <= 186
        return 1.0769e-5 *A - 1.7131e-3
    elseif A > 186 && A <= 208
        return -8.6364e-6 *A + 1.8964e-3
    elseif A > 208
        return 1e-4
    end
end

function Strength_function_old(A)
    if A <= 70
        return 7e-5
    elseif A > 70 && A <= 86
        return 1e-4 *(A *1.875e-2 - 6.125e-1)
    elseif A > 86 && A <= 111
        return 1e-4 
    elseif A > 111 && A <= 121
        return 1e-4 *(-A *2.857e-2 + 4.1714)
    elseif A > 121 && A <= 140
        return 1e-4 
    elseif A > 140 && A <= 144
        return 1e-4 *(A *7.5e-2 - 9.8)
    elseif A > 144
        return 1e-4
    end
end

function slope_intercept(x_1, y_1, x_2, y_2)
    a = (y_1-y_2)/(x_1-x_2)
    b = y_1 - a*x_1
    return a, b
end
slope_intercept(150, 0.9, 167, 0.68)

S0_old = Strength_function_old.(sort(unique(RIPL3_data.A)))
S0_new = Strength_function_new.(sort(unique(RIPL3_data.A)))

pltlog = scatter(
    RIPL3_data.A, RIPL3_data.S_0, 
    yscale=:log10, size = (900, 700), color = :black, xlabel = "A", ylabel = "S₀", 
    framestyle = :box, minorgrid = true, title = "", 
    xlims = (15, 255), ylims = (9e-6, 1e-3),
    xticks = (20:10:250), label="Date RIPL-3", dpi=600
    )
pltlog = plot!(pltlog, sort(unique(RIPL3_data.A)), S0_new, color = :red, label="Parametrizare nouă")
pltlog = plot!(pltlog, sort(unique(RIPL3_data.A)), S0_old, color = :blue, label="Parametrizare veche")

pltlin = scatter(
    RIPL3_data.A, RIPL3_data.S_0, 
    yscale=:identity, size = (900, 700), color = :black, xlabel = "A", ylabel = "S₀", 
    framestyle = :box, minorgrid = true, title = "", 
    xlims = (15, 255), ylims = (0.0, 7e-4),
    xticks = (20:10:250), label="Date RIPL-3", dpi=600, yformatter = :scientific
    )
pltlin = plot!(pltlin, sort(unique(RIPL3_data.A)), S0_new, color = :red, label="Parametrizare nouă")
pltlin = plot!(pltlin, sort(unique(RIPL3_data.A)), S0_old, color = :blue, label="Parametrizare veche")