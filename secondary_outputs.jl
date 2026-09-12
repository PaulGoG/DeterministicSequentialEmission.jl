#=
This part of the program processes the main output so that
it can be compared with experimental data (averaging main output data over experimental yield distributions)
=#
#####
#Compute Y(A,Z,TKE) form experimental Y(A,TKE) data
function Process_yield_data(A_0, fragmdomain::Distribution, dY::DataFrame)
    y_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    #Completes fragmentation domain in case data provided only for HF
    if !isassigned(dY.A[dY.A .< A_0/2], 1)
        for A in sort(unique(dY.A), rev = true)
            for Z in fragmdomain.Z[fragmdomain.A .== A_0 - A]
                P_A_Z = fragmdomain.Value[(fragmdomain.A .== A_0 - A) .& (fragmdomain.Z .== Z)][1]
                for TKE in sort(unique(dY.TKE[dY.A .== A]))
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    push!(y_A_Z_TKE.A, A_0 - A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                end
            end
        end
    end
    for A in sort(unique(dY.A))
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
            for TKE in sort(unique(dY.TKE[dY.A .== A]))
                if !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                end
            end
        end
    end
    #Renormalizes data
    f = 200/sum(y_A_Z_TKE.Value)
    y_A_Z_TKE.Value .*= f
    y_A_Z_TKE.σ .*= f
    return y_A_Z_TKE
end
#Average q(A,Z,TKE) over p(Z,A) so it becomes q(A,TKE)
function Average_over_Z(q_A_Z_TKE, fragmdomain::Distribution)
    q_A_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in sort(unique(q_A_Z_TKE.A))
        for TKE in unique(q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A)])
            Denominator = 0.0
            Numerator = 0.0
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
                    Denominator += P_A_Z
                    Numerator += P_A_Z *value
                end
            end
            if Denominator > 0
                push!(q_A_TKE.A, A)
                push!(q_A_TKE.TKE, TKE)
                push!(q_A_TKE.Value, Numerator/Denominator)
            end
        end
    end
    return q_A_TKE
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(A)
function Average_over_TKE_Z(q_A_Z_TKE, y_A_Z_TKE::Distribution)
    q_A = Distribution_unidym(Int[], Float64[], Float64[])
    for A in sort(unique(q_A_Z_TKE.A))
        Denominator = 0.0
        Numerator = 0.0
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_A.Argument, A)
            push!(q_A.Value, Numerator/Denominator)
        end
    end
    return q_A
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(TKE)
function Average_over_A_Z(q_A_Z_TKE, y_A_Z_TKE::Distribution)
    q_TKE = Distribution_unidym(Float64[], Float64[], Float64[])
    for TKE in sort(unique(q_A_Z_TKE.TKE))
        Denominator = 0.0
        Numerator = 0.0
        for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.TKE .== TKE)])
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_TKE.Argument, TKE)
            push!(q_TKE.Value, Numerator/Denominator)
        end
    end
    return q_TKE
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(Z)
function Average_over_A_TKE(q_A_Z_TKE, y_A_Z_TKE::Distribution)
    q_Z = Distribution_unidym(Float64[], Float64[], Float64[])
    for Z in sort(unique(q_A_Z_TKE.Z))
        Denominator = 0.0
        Numerator = 0.0
        for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.Z .== Z)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_Z.Argument, Z)
            push!(q_Z.Value, Numerator/Denominator)
        end
    end
    return q_Z
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) to get average value <q>
function Average_value(q_A_Z_TKE, y_A_Z_TKE::Distribution, mass_number_range)
    Denominator = 0.0
    Numerator = 0.0
    for A in mass_number_range
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in unique(q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)])
                value = sum(q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)])
                if !isnan(value)
                    n = length(q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)])
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE *n
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
    end
    return Numerator/Denominator
end
#Average q(Argument) over Y(Argument) to get average value <q>
function Average_value(q, y::Distribution_unidym, Argument_range)
    Denominator = 0.0
    Numerator = 0.0
    for Argument in Argument_range
        if isassigned(q.Value[q.Argument .== Argument], 1) && isassigned(y.Value[y.Argument .== Argument], 1)
            value = q.Value[q.Argument .== Argument][1]
            Y = y.Value[y.Argument .== Argument][1]
            Denominator += Y
            Numerator += Y *value
        end
    end
    return Numerator/Denominator
end
#Get L-H pair value from given q(A) distribution and A_H 
function Pair_value(q_A::Distribution_unidym, A_0, A_H)
    val_L = q_A.Value[q_A.Argument .== A_0 - A_H][1]
    val_H = q_A.Value[q_A.Argument .== A_H][1]
    return val_L + val_H
end
#Obtain Yp(Z,Aₚ,TKE), Yp(Z,Aₚ) distributions from Y(A,Z,TKE) & n(A,Z,TKE)
function Yield_post_neutron(y_A_Z_TKE::Distribution, n_A_Z_TKE)
    yₚ_Aₚ_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    yₚ_Aₚ_Z = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(y_A_Z_TKE.A)
        for Z in unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A])
            for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z)]
                Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                σY_A_Z_TKE = y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                if isassigned(n_A_Z_TKE.Value[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.Z .== Z) .& (n_A_Z_TKE.TKE .== TKE)], 1)
                    n = n_A_Z_TKE.Value[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.Z .== Z) .& (n_A_Z_TKE.TKE .== TKE)][1]
                    Aₚ = A - n
                    if !isassigned(yₚ_Aₚ_Z_TKE.Value[(yₚ_Aₚ_Z_TKE.A .== Aₚ) .& (yₚ_Aₚ_Z_TKE.Z .== Z) .& (yₚ_Aₚ_Z_TKE.TKE .== TKE)], 1)
                        push!(yₚ_Aₚ_Z_TKE.A, Aₚ)
                        push!(yₚ_Aₚ_Z_TKE.Z, Z)
                        push!(yₚ_Aₚ_Z_TKE.TKE, TKE)
                        push!(yₚ_Aₚ_Z_TKE.Value, Y_A_Z_TKE)
                        push!(yₚ_Aₚ_Z_TKE.σ, σY_A_Z_TKE)
                    else
                        yₚ_Aₚ_Z_TKE.Value[(yₚ_Aₚ_Z_TKE.A .== Aₚ) .& (yₚ_Aₚ_Z_TKE.Z .== Z) .& (yₚ_Aₚ_Z_TKE.TKE .== TKE)] .+= Y_A_Z_TKE
                        yₚ_Aₚ_Z_TKE.σ[(yₚ_Aₚ_Z_TKE.A .== Aₚ) .& (yₚ_Aₚ_Z_TKE.Z .== Z) .& (yₚ_Aₚ_Z_TKE.TKE .== TKE)] .= sqrt(sum(yₚ_Aₚ_Z_TKE.σ[(yₚ_Aₚ_Z_TKE.A .== Aₚ) .& (yₚ_Aₚ_Z_TKE.Z .== Z) .& (yₚ_Aₚ_Z_TKE.TKE .== TKE)].^2) + σY_A_Z_TKE^2)
                    end
                end
            end
        end
    end
    f = 200/sum(yₚ_Aₚ_Z_TKE.Value)
    yₚ_Aₚ_Z_TKE.Value .*= f
    yₚ_Aₚ_Z_TKE.σ .*= f
    for A in sort(unique(yₚ_Aₚ_Z_TKE.A))
        for Z in sort(unique(yₚ_Aₚ_Z_TKE.Z[(yₚ_Aₚ_Z_TKE.A .== A)]))
            Yₚ_Aₚ_Z = sum(yₚ_Aₚ_Z_TKE.Value[(yₚ_Aₚ_Z_TKE.A .== A) .& (yₚ_Aₚ_Z_TKE.Z .== Z)])
            σYₚ_Aₚ_Z = sqrt(sum(yₚ_Aₚ_Z_TKE.σ[(yₚ_Aₚ_Z_TKE.A .== A) .& (yₚ_Aₚ_Z_TKE.Z .== Z)].^2))
            push!(yₚ_Aₚ_Z.A, A)
            push!(yₚ_Aₚ_Z.Z, Z)
            push!(yₚ_Aₚ_Z.Value, Yₚ_Aₚ_Z)
            push!(yₚ_Aₚ_Z.σ, σYₚ_Aₚ_Z)
        end
    end
    return yₚ_Aₚ_Z_TKE, yₚ_Aₚ_Z 
end
#Compute KEₚ(A,Z,TKE) and TKEₚ(A_H,Z_H,TKE)
function Kinetic_Energy_post_neutron(A_0, Z_0, A_H_range, n_A_Z_TKE)
    keₚ_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    tkeₚ_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(n_A_Z_TKE.A)
        for TKE in unique(n_A_Z_TKE.TKE[n_A_Z_TKE.A .== A])
            KE = TKE *(A_0 - A)/A_0
            for Z in n_A_Z_TKE.Z[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.TKE .== TKE)]
                n = n_A_Z_TKE.Value[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.Z .== Z) .& (n_A_Z_TKE.TKE .== TKE)][1]
                KEp = KE *(A - n)/A
                push!(keₚ_A_Z_TKE.A, A)
                push!(keₚ_A_Z_TKE.Z, Z)
                push!(keₚ_A_Z_TKE.TKE, TKE)
                push!(keₚ_A_Z_TKE.Value, KEp)
            end
        end
    end
    for A in A_H_range
        for Z in unique(keₚ_A_Z_TKE.Z[keₚ_A_Z_TKE.A .== A])
            for TKE in unique(keₚ_A_Z_TKE.TKE[(keₚ_A_Z_TKE.A .== A) .& (keₚ_A_Z_TKE.Z .== Z)])
                if !isassigned(tkeₚ_A_Z_TKE.Value[(tkeₚ_A_Z_TKE.A .== A) .& (tkeₚ_A_Z_TKE.Z .== Z) .& (tkeₚ_A_Z_TKE.TKE .== TKE)], 1)
                    KEp_L = keₚ_A_Z_TKE.Value[(keₚ_A_Z_TKE.A .== A_0 - A) .& (keₚ_A_Z_TKE.Z .== Z_0 - Z) .& (keₚ_A_Z_TKE.TKE .== TKE)][1]
                    KEp_H = keₚ_A_Z_TKE.Value[(keₚ_A_Z_TKE.A .== A) .& (keₚ_A_Z_TKE.Z .== Z) .& (keₚ_A_Z_TKE.TKE .== TKE)][1] 
                    TKEp = KEp_L + KEp_H
                    push!(tkeₚ_A_Z_TKE.A, A)
                    push!(tkeₚ_A_Z_TKE.Z, Z)
                    push!(tkeₚ_A_Z_TKE.TKE, TKE)
                    push!(tkeₚ_A_Z_TKE.Value, TKEp)
                end
            end
        end
    end
    return keₚ_A_Z_TKE, tkeₚ_A_Z_TKE
end
function Average_Kinetic_Energy_post_neutron(A_H_min, keₚ_A_Z_TKE, tkeₚ_A_Z_TKE, y_A_Z_TKE, n_A_Z_TKE)
    keₚ_Ap = Distribution_unidym(Int[], Float64[], Float64[])
    tkeₚ_AHp = Distribution_unidym(Int[], Float64[], Float64[])
    for A in unique(y_A_Z_TKE.A)
        for Z in unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A])
            for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z)]
                if isassigned(keₚ_A_Z_TKE.Value[(keₚ_A_Z_TKE.A .== A) .& (keₚ_A_Z_TKE.Z .== Z) .& (keₚ_A_Z_TKE.TKE .== TKE)], 1)
                    KEp = keₚ_A_Z_TKE.Value[(keₚ_A_Z_TKE.A .== A) .& (keₚ_A_Z_TKE.Z .== Z) .& (keₚ_A_Z_TKE.TKE .== TKE)][1]
                    Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                    n = n_A_Z_TKE.Value[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.Z .== Z) .& (n_A_Z_TKE.TKE .== TKE)][1]
                    Aₚ = A - n
                    if !isassigned(keₚ_Ap.Value[keₚ_Ap.Argument .== Aₚ], 1)
                        push!(keₚ_Ap.Argument, Aₚ)
                        push!(keₚ_Ap.Value, KEp *Y_A_Z_TKE)
                        push!(keₚ_Ap.σ, Y_A_Z_TKE)
                    else
                        keₚ_Ap.Value[keₚ_Ap.Argument .== Aₚ] .+= KEp *Y_A_Z_TKE
                        keₚ_Ap.σ[keₚ_Ap.Argument .== Aₚ] .+= Y_A_Z_TKE
                    end
                    if A >= A_H_min
                        TKEp = tkeₚ_A_Z_TKE.Value[(tkeₚ_A_Z_TKE.A .== A) .& (tkeₚ_A_Z_TKE.Z .== Z) .& (tkeₚ_A_Z_TKE.TKE .== TKE)][1]
                        if !isassigned(tkeₚ_AHp.Value[tkeₚ_AHp.Argument .== Aₚ], 1)
                            push!(tkeₚ_AHp.Argument, Aₚ)
                            push!(tkeₚ_AHp.Value, TKEp *Y_A_Z_TKE)
                            push!(tkeₚ_AHp.σ, Y_A_Z_TKE)
                        else
                            tkeₚ_AHp.Value[tkeₚ_AHp.Argument .== Aₚ] .+= TKEp *Y_A_Z_TKE
                            tkeₚ_AHp.σ[tkeₚ_AHp.Argument .== Aₚ] .+= Y_A_Z_TKE
                        end
                    end
                end
            end
        end
    end
    keₚ_Ap.Value ./= keₚ_Ap.σ
    tkeₚ_AHp.Value ./= tkeₚ_AHp.σ
    empty!(keₚ_Ap.σ)
    empty!(tkeₚ_AHp.σ)
    Sort_q_Argument(keₚ_Ap)
    Sort_q_Argument(tkeₚ_AHp)
    return keₚ_Ap, tkeₚ_AHp
end
#Compute RT(A_H,Z_H,TKE)
function Ratio_of_Temperatures(A_0, Z_0, A_H_range, fragmdomain, E_exi_A_Z_TKE::Distribution, density_parameter_type, density_parameter_data)
    T0_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    RT_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in A_H_range
        for Z in unique(E_exi_A_Z_TKE.Z[E_exi_A_Z_TKE.A .== A])
            a_L = density_parameter(density_parameter_type, A_0 - A, Z_0 - Z, density_parameter_data)
            a_H = density_parameter(density_parameter_type, A, Z, density_parameter_data) 
            for TKE in E_exi_A_Z_TKE.TKE[(E_exi_A_Z_TKE.A .== A) .& (E_exi_A_Z_TKE.Z .== Z)]
                if !isassigned(RT_A_Z_TKE.Value[(RT_A_Z_TKE.A .== A) .& (RT_A_Z_TKE.Z .== Z) .& (RT_A_Z_TKE.TKE .== TKE)], 1)
                    E_exi_L = E_exi_A_Z_TKE.Value[(E_exi_A_Z_TKE.A .== A_0 - A) .& (E_exi_A_Z_TKE.Z .== Z_0 - Z) .& (E_exi_A_Z_TKE.TKE .== TKE)][1]
                    E_exi_H = E_exi_A_Z_TKE.Value[(E_exi_A_Z_TKE.A .== A) .& (E_exi_A_Z_TKE.Z .== Z) .& (E_exi_A_Z_TKE.TKE .== TKE)][1] 
                    T_L = sqrt(E_exi_L/a_L)
                    T_H = sqrt(E_exi_H/a_H)
                    RT = T_L/T_H

                    push!(RT_A_Z_TKE.A, A)
                    push!(RT_A_Z_TKE.Z, Z)
                    push!(RT_A_Z_TKE.TKE, TKE)
                    push!(RT_A_Z_TKE.Value, RT)

                    push!(T0_A_Z_TKE.A, A)
                    push!(T0_A_Z_TKE.Z, Z)
                    push!(T0_A_Z_TKE.TKE, TKE)
                    push!(T0_A_Z_TKE.Value, T_H)
                    push!(T0_A_Z_TKE.A, A_0 - A)
                    push!(T0_A_Z_TKE.Z, Z_0 - Z)
                    push!(T0_A_Z_TKE.TKE, TKE)
                    push!(T0_A_Z_TKE.Value, T_L)
                end
            end
        end
    end
    Sort_q_A_Z_TKE(T0_A_Z_TKE, fragmdomain)
    return RT_A_Z_TKE, T0_A_Z_TKE
end
#Obtain vectorized singular distributions Y(A), Y(N), Y(Z), Y(TKE), TKE(AH) & KE(A) from Y(A,Z,TKE)
function Singular_yield_distributions(y_A_Z_TKE::Distribution, A_0, A_H_min)
    y_A = Distribution_unidym(Int[], Float64[], Float64[])
    y_Z = Distribution_unidym(Int[], Float64[], Float64[])
    y_N = Distribution_unidym(Int[], Float64[], Float64[])
    y_TKE = Distribution_unidym(Float64[], Float64[], Float64[])
    tke_AH = Distribution_unidym(Int[], Float64[], Float64[])
    ke_A = Distribution_unidym(Int[], Float64[], Float64[])
    #Y(A)
    for A in sort(unique(y_A_Z_TKE.A))
        Y_A = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .== A])
        σY_A = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.A .== A] .^2))
        push!(y_A.Argument, A)
        push!(y_A.Value, Y_A)
        push!(y_A.σ, σY_A)
    end
    #Y(Z)
    for Z in sort(unique(y_A_Z_TKE.Z))
        Y_Z = sum(y_A_Z_TKE.Value[y_A_Z_TKE.Z .== Z])
        σY_Z = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.Z .== Z] .^2))
        push!(y_Z.Argument, Z)
        push!(y_Z.Value, Y_Z)
        push!(y_Z.σ, σY_Z)
    end
    #Y(N)
    for A in sort(unique(y_A_Z_TKE.A))
        for Z in unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A])
            N = A - Z
            if !isassigned(y_N.Value[y_N.Argument .== N], 1)
                Y_N = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .- y_A_Z_TKE.Z .== N])
                σY_N = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.A .- y_A_Z_TKE.Z .== N].^2))
                push!(y_N.Argument, N)
                push!(y_N.Value, Y_N)
                push!(y_N.σ, σY_N)
            end
        end
    end
    #Y(TKE)
    for TKE in sort(unique(y_A_Z_TKE.TKE))
        Y_TKE = sum(y_A_Z_TKE.Value[y_A_Z_TKE.TKE .== TKE])
        σY_TKE = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.TKE .== TKE] .^2))
        push!(y_TKE.Argument, TKE)
        push!(y_TKE.Value, Y_TKE)
        push!(y_TKE.σ, σY_TKE)
    end
    #TKE(AH)
    for A_H in sort(unique(y_A_Z_TKE.A[y_A_Z_TKE.A .>= A_H_min]))
        Numerator = 0.0
        𝚺_σ² = 0.0
        Denominator = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .== A_H])
        for TKE in sort(unique(y_A_Z_TKE.TKE[y_A_Z_TKE.A .== A_H]))
            Y_A_TKE = sum(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.TKE .== TKE)])
            Numerator += TKE *Y_A_TKE
        end
        TKE_AH = Numerator/Denominator
        for TKE in sort(unique(y_A_Z_TKE.TKE[y_A_Z_TKE.A .== A_H]))
            σY_A_TKE = sqrt(sum(y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.TKE .== TKE)] .^2))
            𝚺_σ² += (TKE - TKE_AH)^2 *σY_A_TKE^2
        end
        push!(tke_AH.Argument, A_H)
        push!(tke_AH.Value, TKE_AH)
        push!(tke_AH.σ, sqrt(𝚺_σ²)/Denominator)
    end  
    #KE(A)
    for A_H in sort(unique(tke_AH.Argument), rev = true)
        TKE_A = tke_AH.Value[tke_AH.Argument .== A_H][1]
        σTKE_A = tke_AH.σ[tke_AH.Argument .== A_H][1]
        KE_A = TKE_A *A_H /A_0
        push!(ke_A.Argument, A_0 - A_H)
        push!(ke_A.Value, KE_A)
        push!(ke_A.σ, KE_A *σTKE_A /TKE_A)
    end
    for A_H in unique(tke_AH.Argument)
        if A_H != A_0 - A_H
            TKE_A = tke_AH.Value[tke_AH.Argument .== A_H][1]
            σTKE_A = tke_AH.σ[tke_AH.Argument .== A_H][1]
            KE_A = TKE_A *(A_0 - A_H) /A_0
            push!(ke_A.Argument, A_H)
            push!(ke_A.Value, KE_A)
            push!(ke_A.σ, KE_A *σTKE_A /TKE_A)
        end
    end
    Sort_q_Argument(y_A)
    Sort_q_Argument(y_Z)
    Sort_q_Argument(y_N)
    Sort_q_Argument(y_TKE)
    Sort_q_Argument(tke_AH)
    Sort_q_Argument(ke_A)
    return y_A, y_Z, y_N, y_TKE, tke_AH, ke_A
end
#Get <Argument> average value from Yield(Argument) singular yield distribution
function Average_yield_argument(yield::Distribution_unidym, argument_range)
    Numerator = 0.0
    Denominator = 0.0
    𝚺_σ² = 0.0
    for Argument in argument_range
        Y = yield.Value[yield.Argument .== Argument][1]
        Numerator += Argument *Y
        Denominator += Y
    end
    Average_arg = Numerator/Denominator
    for Argument in argument_range
        σY = yield.σ[yield.Argument .== Argument][1]
        𝚺_σ² += (Argument - Average_arg)^2 *σY^2
    end    
    return Average_arg, sqrt(𝚺_σ²)/Denominator
end
#Obtain vectorized distributions Q(A,Z), Q(AH), TXE(A,Z,TKE)
function Vectorized_TXE_Q(A_0, Z_0, fission_type::String, E_incident, y_A_Z_TKE::Distribution, A_H_range, dm::DataFrame)
    q_A_Z =  Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Q_AH = Distribution_unidym(Int[], Float64[], Float64[])
    txe_A_Z_TKE =  Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    E_CN = Compound_nucleus_energy(fission_type, A_0, Z_0, E_incident, dm)
    for A_H in A_H_range
        Denominator_Q_AH = 0.0
        Numerator_Q_AH = 0.0
        for Z_H in unique(y_A_Z_TKE.Z[(y_A_Z_TKE.A .== A_H)])
            Q_A_Z = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
            if !isnan(Q_A_Z[1])
                push!(q_A_Z.A, A_H)
                push!(q_A_Z.Z, Z_H)
                push!(q_A_Z.Value, Q_A_Z[1])
                push!(q_A_Z.σ, Q_A_Z[2])
                for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)]
                    TXE_A_Z_TKE = Total_excitation_energy(Q_A_Z[1], Q_A_Z[2], TKE, 0.0, E_CN[1], E_CN[2])
                    if !isnan(TXE_A_Z_TKE[1])
                        push!(txe_A_Z_TKE.A, A_H)
                        push!(txe_A_Z_TKE.Z, Z_H)
                        push!(txe_A_Z_TKE.TKE, TKE)
                        push!(txe_A_Z_TKE.Value, TXE_A_Z_TKE[1])
                        push!(txe_A_Z_TKE.σ, TXE_A_Z_TKE[2])
                    end
                end
                Y_A_Z = sum(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)])
                Numerator_Q_AH += Y_A_Z *Q_A_Z[1]
                Denominator_Q_AH += Y_A_Z
            end
        end
        if Denominator_Q_AH > 0
            push!(Q_AH.Argument, A_H)
            push!(Q_AH.Value, Numerator_Q_AH/Denominator_Q_AH)
        end
    end
    return txe_A_Z_TKE, q_A_Z, Q_AH
end
#Compute probabilities of occurance P(q) of any vectorized generic quantity q(A,Z,TKE) with value step of Δq
function Probability_of_occurrence(q_A_Z_TKE, y_A_Z_TKE::Distribution, Δq::Number)
    P = Distribution_unidym(Float64[], Float64[], Float64[])
    q_values = filter(!isnan, q_A_Z_TKE.Value)
    if isassigned(q_values, 1)
        lower_q = floor(minimum(q_values))
        upper_q = lower_q + Δq
        while lower_q <= maximum(q_values)
            Frequency = 0.0
            for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.Value .>= lower_q) .& (q_A_Z_TKE.Value .< upper_q)])
                for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Value .>= lower_q) .& (q_A_Z_TKE.Value .< upper_q)])
                    for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z) .& (q_A_Z_TKE.Value .>= lower_q) .& (q_A_Z_TKE.Value .< upper_q)]
                        if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                            Frequency += y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        end
                    end
                end
            end
            push!(P.Argument, lower_q)
            push!(P.Value, Frequency)
            lower_q = upper_q
            upper_q += Δq
        end
        #f = 100/sum(P.Value)
        #P.Value .*= f
    else
        push!(P.Argument, NaN)
        push!(P.Value, NaN)
    end
    return P
end
#Sort q(Argument) by ascending argument for plotting
function Sort_q_Argument(q)
    aux_Argument = sort(unique(q.Argument))
    aux_Value = [first(q.Value[q.Argument .== Argument]) for Argument in aux_Argument]
    for index in eachindex(aux_Argument)
        q.Argument[index] = aux_Argument[index]
        q.Value[index] = aux_Value[index]
    end
    if isassigned(filter(!isnan, q.σ), 1)
        aux_σ = [first(q.σ[q.Argument .== Argument]) for Argument in aux_Argument]
        for index in eachindex(aux_Argument)
            q.σ[index] = aux_σ[index]
        end
    end
end
#####
y_A_Z_TKE = Process_yield_data(A₀, fragmdomain, Yield_data)

Write_seq_output_Y(A₀, Z₀, A_H_min, A_H_max, No_ZperA, tkerange, 
y_A_Z_TKE, E_excitation, Raw_output_datafile, density_parameter_type, density_parameter_data, 
fissionant_nucleus_identifier, file_output_identifier, mass_excess_filename, 
txe_partitioning_type, txe_partitioning_data, evaporation_cs_type, dmass_excess)
println("*averaging data over $yield_distribution_filename provided Yield distribution") 

if secondary_output_Yield
    y_A, y_Z, y_N, y_TKE, tke_AH, ke_A = Singular_yield_distributions(y_A_Z_TKE, A₀, A_H_min)
    if !isdir("$(file_output_identifier)_output_data/Yield/")
        mkdir("$(file_output_identifier)_output_data/Yield/")
    end
    if isassigned(filter(!isnan, y_A_Z_TKE.σ), 1)
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_A_Z_TKE_$(file_output_identifier).dat", 
            DataFrame(A = y_A_Z_TKE.A, Z = y_A_Z_TKE.Z, TKE = y_A_Z_TKE.TKE, Y = y_A_Z_TKE.Value, σ = y_A_Z_TKE.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_A_$(file_output_identifier).dat", 
            DataFrame(A = y_A.Argument, Y = y_A.Value, σ = y_A.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_Z_$(file_output_identifier).dat", 
            DataFrame(Z = y_Z.Argument, Y = y_Z.Value, σ = y_Z.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_N_$(file_output_identifier).dat", 
            DataFrame(N = y_N.Argument, Y = y_N.Value, σ = y_N.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_TKE_$(file_output_identifier).dat", 
            DataFrame(TKE = y_TKE.Argument, Y = y_TKE.Value, σ = y_TKE.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_TKE_AH_$(file_output_identifier).dat", 
            DataFrame(A_H = tke_AH.Argument, TKE = tke_AH.Value, σ = tke_AH.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_KE_A_$(file_output_identifier).dat", 
            DataFrame(A = ke_A.Argument, KE = ke_A.Value, σ = ke_A.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
    else
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_A_Z_TKE_$(file_output_identifier).dat", 
            DataFrame(A = y_A_Z_TKE.A, Z = y_A_Z_TKE.Z, TKE = y_A_Z_TKE.TKE, Y = y_A_Z_TKE.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_A_$(file_output_identifier).dat", 
            DataFrame(A = y_A.Argument, Y = y_A.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_Z_$(file_output_identifier).dat", 
            DataFrame(Z = y_Z.Argument, Y = y_Z.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_N_$(file_output_identifier).dat", 
            DataFrame(N = y_N.Argument, Y = y_N.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_Y_TKE_$(file_output_identifier).dat", 
            DataFrame(TKE = y_TKE.Argument, Y = y_TKE.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_TKE_AH_$(file_output_identifier).dat", 
            DataFrame(A_H = tke_AH.Argument, TKE = tke_AH.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Yield/$(fissionant_nucleus_identifier)_KE_A_$(file_output_identifier).dat", 
            DataFrame(A = ke_A.Argument, KE = ke_A.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
    end
end
if secondary_output_nu
    if !isdir("$(file_output_identifier)_output_data/nu/")
        mkdir("$(file_output_identifier)_output_data/nu/")
    end
    ν_A_Z_TKE = Neutron_multiplicity_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence
        )
    )
    ν_Pair_A_Z_TKE = Neutron_multiplicity_Pair_A_Z_TKE(A₀, Z₀, A_H_range, ν_A_Z_TKE)
    max_seq_A_Z_TKE = Maximum_sequences_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence
        )
    )
    ν_A_TKE = Average_over_Z(ν_A_Z_TKE, fragmdomain)
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_nu_A_TKE_$(file_output_identifier).dat", 
        DataFrame(A = ν_A_TKE.A, TKE = ν_A_TKE.TKE, nu = ν_A_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_A = Average_over_TKE_Z(ν_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_nu_A_$(file_output_identifier).dat", 
        DataFrame(A = ν_A.Argument, nu = ν_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_AH_Pair = DataFrame(
        Argument =  ν_A.Argument[ν_A.Argument .>= A_H_min],
        Value = [Pair_value(ν_A, A₀, A_H) for A_H in ν_A.Argument[ν_A.Argument .>= A_H_min]]
    )
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_nu_AH_Pair_$(file_output_identifier).dat", 
        DataFrame(A = ν_AH_Pair.Argument, nu_Pair = ν_AH_Pair.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_L_TKE = Average_over_A_Z(
        DataFrame(
        A = ν_A_Z_TKE.A[ν_A_Z_TKE.A .<= A_H_min], 
        Z = ν_A_Z_TKE.Z[ν_A_Z_TKE.A .<= A_H_min], 
        TKE = ν_A_Z_TKE.TKE[ν_A_Z_TKE.A .<= A_H_min],
        Value = ν_A_Z_TKE.Value[ν_A_Z_TKE.A .<= A_H_min]
    ), y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_nu_LF_TKE_$(file_output_identifier).dat", 
        DataFrame(TKE = ν_L_TKE.Argument, nu = ν_L_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_H_TKE = Average_over_A_Z(
        DataFrame(
        A = ν_A_Z_TKE.A[ν_A_Z_TKE.A .>= A_H_min], 
        Z = ν_A_Z_TKE.Z[ν_A_Z_TKE.A .>= A_H_min], 
        TKE = ν_A_Z_TKE.TKE[ν_A_Z_TKE.A .>= A_H_min],
        Value = ν_A_Z_TKE.Value[ν_A_Z_TKE.A .>= A_H_min]
    ), y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_nu_HF_TKE_$(file_output_identifier).dat", 
        DataFrame(TKE = ν_H_TKE.Argument, nu = ν_H_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_Pair_TKE = Average_over_A_Z(ν_Pair_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_nu_Pair_TKE_$(file_output_identifier).dat", 
        DataFrame(TKE = ν_Pair_TKE.Argument, nu = ν_Pair_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_ν_Pair = Probability_of_occurrence(ν_Pair_A_Z_TKE, y_A_Z_TKE, 1)
    CSV.write(
        "$(file_output_identifier)_output_data/nu/$(fissionant_nucleus_identifier)_P_nu_Pair_$(file_output_identifier).dat", 
        DataFrame(nu = probability_ν_Pair.Argument, P = probability_ν_Pair.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    if secondary_output_Ap
        if !isdir("$(file_output_identifier)_output_data/Post_neutron/")
            mkdir("$(file_output_identifier)_output_data/Post_neutron/")
        end
        if !isdir("$(file_output_identifier)_output_data/Post_neutron/Yield_Ap_Z/")
            mkdir("$(file_output_identifier)_output_data/Post_neutron/Yield_Ap_Z/")
        end
        yp_Ap_Z_TKE, yp_Ap_Z = Yield_post_neutron(y_A_Z_TKE, max_seq_A_Z_TKE)
        yp_Ap, placeholder1, yp_Np, placeholder2, placeholder3, placeholder4 = Singular_yield_distributions(yp_Ap_Z_TKE, A₀, A_H_min)
        if isassigned(filter(!isnan, yp_Ap_Z_TKE.σ), 1)
            CSV.write(
                "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_Yp_Ap_$(file_output_identifier).dat", 
                DataFrame(Ap = yp_Ap.Argument, Yp = yp_Ap.Value, σ = yp_Ap.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_Yp_Ap_Z_$(file_output_identifier).dat", 
                DataFrame(Ap = yp_Ap_Z.A, Z = yp_Ap_Z.Z, Yp = yp_Ap_Z.Value, σ = yp_Ap_Z.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_Yp_Np_$(file_output_identifier).dat", 
                DataFrame(Np = yp_Np.Argument, Yp = yp_Np.Value, σ = yp_Np.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            for Z in unique(yp_Ap_Z.Z)
                CSV.write(
                    "$(file_output_identifier)_output_data/Post_neutron/Yield_Ap_Z/$(fissionant_nucleus_identifier)_Yp_Ap_$(Z)_$(file_output_identifier).dat", 
                    DataFrame(Ap = yp_Ap_Z.A[yp_Ap_Z.Z .== Z], Yp = yp_Ap_Z.Value[yp_Ap_Z.Z .== Z], σ = yp_Ap_Z.σ[yp_Ap_Z.Z .== Z]), 
                    writeheader=true, newline="\r\n", delim=' '
                )
            end
        else
            CSV.write(
                "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_Yp_Ap_$(file_output_identifier).dat", 
                DataFrame(Ap = yp_Ap.Argument, Yp = yp_Ap.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_Yp_Ap_Z_$(file_output_identifier).dat", 
                DataFrame(Ap = yp_Ap_Z.A, Z = yp_Ap_Z.Z, Yp = yp_Ap_Z.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_Yp_Np_$(file_output_identifier).dat", 
                DataFrame(Np = yp_Np.Argument, Yp = yp_Np.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            for Z in unique(yp_Ap_Z.Z)
                CSV.write(
                    "$(file_output_identifier)_output_data/Post_neutron/Yield_Ap_Z/$(fissionant_nucleus_identifier)_Yp_Ap_$(Z)_$(file_output_identifier).dat", 
                    DataFrame(Ap = yp_Ap_Z.A[yp_Ap_Z.Z .== Z], Yp = yp_Ap_Z.Value[yp_Ap_Z.Z .== Z]), 
                    writeheader=true, newline="\r\n", delim=' '
                )
            end
        end
        kep_A_Z_TKE, tkep_A_Z_TKE = Kinetic_Energy_post_neutron(A₀, Z₀, A_H_range, max_seq_A_Z_TKE)
        kep_Z = Average_over_A_TKE(kep_A_Z_TKE, y_A_Z_TKE)
        kep_Ap, tkep_AHp = Average_Kinetic_Energy_post_neutron(A_H_min, kep_A_Z_TKE, tkep_A_Z_TKE, y_A_Z_TKE, max_seq_A_Z_TKE)
        CSV.write(
            "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_TKEp_AHp_$(file_output_identifier).dat", 
            DataFrame(Ap = tkep_AHp.Argument, TKEp = tkep_AHp.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_KEp_Ap_$(file_output_identifier).dat", 
            DataFrame(Ap = kep_Ap.Argument, KEp = kep_Ap.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "$(file_output_identifier)_output_data/Post_neutron/$(fissionant_nucleus_identifier)_KEp_Z_$(file_output_identifier).dat", 
            DataFrame(Z = kep_Z.Argument, KEp = kep_Z.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
    end
    if secondary_output_Tₖ
        if !isdir("$(file_output_identifier)_output_data/P_T_k/")
            mkdir("$(file_output_identifier)_output_data/P_T_k/")
        end
        if !isdir("$(file_output_identifier)_output_data/T_k_A/")
            mkdir("$(file_output_identifier)_output_data/T_k_A/")
        end
        for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
            Tₖ_A_Z_TKE = DataFrame(
                A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Value = Raw_output_datafile.Tₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
            )
            Tₖ_A = Average_over_TKE_Z(Tₖ_A_Z_TKE, y_A_Z_TKE)
            CSV.write(
                "$(file_output_identifier)_output_data/T_k_A/$(fissionant_nucleus_identifier)_T_$(k)_A_$(file_output_identifier).dat", 
                DataFrame(A = Tₖ_A.Argument, T = Tₖ_A.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_Tₖ = Probability_of_occurrence(Tₖ_A_Z_TKE, y_A_Z_TKE, ΔTₖ)
            CSV.write(
                "$(file_output_identifier)_output_data/P_T_k/$(fissionant_nucleus_identifier)_P_T_$(k)_$(file_output_identifier).dat", 
                DataFrame(T = probability_Tₖ.Argument, P = probability_Tₖ.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_Tₖ_L = Probability_of_occurrence(
                DataFrame(
                    A = Tₖ_A_Z_TKE.A[Tₖ_A_Z_TKE.A .<= A_H_min], 
                    Z = Tₖ_A_Z_TKE.Z[Tₖ_A_Z_TKE.A .<= A_H_min], 
                    TKE = Tₖ_A_Z_TKE.TKE[Tₖ_A_Z_TKE.A .<= A_H_min],
                    Value = Tₖ_A_Z_TKE.Value[Tₖ_A_Z_TKE.A .<= A_H_min]
                ), y_A_Z_TKE, ΔTₖ
                )
            CSV.write(
                "$(file_output_identifier)_output_data/P_T_k/$(fissionant_nucleus_identifier)_P_T_$(k)_LF_$(file_output_identifier).dat", 
                DataFrame(T = probability_Tₖ_L.Argument, P = probability_Tₖ_L.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_Tₖ_H = Probability_of_occurrence(
                DataFrame(
                    A = Tₖ_A_Z_TKE.A[Tₖ_A_Z_TKE.A .>= A_H_min], 
                    Z = Tₖ_A_Z_TKE.Z[Tₖ_A_Z_TKE.A .>= A_H_min], 
                    TKE = Tₖ_A_Z_TKE.TKE[Tₖ_A_Z_TKE.A .>= A_H_min],
                    Value = Tₖ_A_Z_TKE.Value[Tₖ_A_Z_TKE.A .>= A_H_min]
                ), y_A_Z_TKE, ΔTₖ
                )
            CSV.write(
                "$(file_output_identifier)_output_data/P_T_k/$(fissionant_nucleus_identifier)_P_T_$(k)_HF_$(file_output_identifier).dat", 
                DataFrame(T = probability_Tₖ_H.Argument, P = probability_Tₖ_H.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            if secondary_output_Eᵣ
                if !isdir("$(file_output_identifier)_output_data/P_Er/")
                    mkdir("$(file_output_identifier)_output_data/P_Er/")
                end
                aₖ = copy(Raw_output_datafile.aₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)])
                Eᵣ_A_Z_TKE = copy(Tₖ_A_Z_TKE)
                Eᵣ_A_Z_TKE.Value .= Energy_FermiGas.(aₖ, Eᵣ_A_Z_TKE.Value)
                probability_Eᵣ = Probability_of_occurrence(Eᵣ_A_Z_TKE, y_A_Z_TKE, ΔEᵣ)
                CSV.write(
                    "$(file_output_identifier)_output_data/P_Er/$(fissionant_nucleus_identifier)_P_Er_$(k)_$(file_output_identifier).dat", 
                    DataFrame(Eᵣ = probability_Eᵣ.Argument, P = probability_Eᵣ.Value), 
                    writeheader=true, newline="\r\n", delim=' '
                )
                probability_Eᵣ_L = Probability_of_occurrence(
                    DataFrame(
                        A = Eᵣ_A_Z_TKE.A[Eᵣ_A_Z_TKE.A .<= A_H_min], 
                        Z = Eᵣ_A_Z_TKE.Z[Eᵣ_A_Z_TKE.A .<= A_H_min], 
                        TKE = Eᵣ_A_Z_TKE.TKE[Eᵣ_A_Z_TKE.A .<= A_H_min],
                        Value = Eᵣ_A_Z_TKE.Value[Eᵣ_A_Z_TKE.A .<= A_H_min]
                    ), y_A_Z_TKE, ΔEᵣ
                    )
                CSV.write(
                    "$(file_output_identifier)_output_data/P_Er/$(fissionant_nucleus_identifier)_P_Er_$(k)_LF_$(file_output_identifier).dat", 
                    DataFrame(Eᵣ = probability_Eᵣ_L.Argument, P = probability_Eᵣ_L.Value), 
                    writeheader=true, newline="\r\n", delim=' '
                )
                probability_Eᵣ_H = Probability_of_occurrence(
                    DataFrame(
                        A = Eᵣ_A_Z_TKE.A[Eᵣ_A_Z_TKE.A .>= A_H_min], 
                        Z = Eᵣ_A_Z_TKE.Z[Eᵣ_A_Z_TKE.A .>= A_H_min], 
                        TKE = Eᵣ_A_Z_TKE.TKE[Eᵣ_A_Z_TKE.A .>= A_H_min],
                        Value = Eᵣ_A_Z_TKE.Value[Eᵣ_A_Z_TKE.A .>= A_H_min]
                    ), y_A_Z_TKE, ΔEᵣ
                    )
                CSV.write(
                    "$(file_output_identifier)_output_data/P_Er/$(fissionant_nucleus_identifier)_P_Er_$(k)_HF_$(file_output_identifier).dat", 
                    DataFrame(Eᵣ = probability_Eᵣ_H.Argument, P = probability_Eᵣ_H.Value), 
                    writeheader=true, newline="\r\n", delim=' '
                )
            end      
        end
    end
    if secondary_output_avg_εₖ
        if !isdir("$(file_output_identifier)_output_data/P_avgE_k/")
            mkdir("$(file_output_identifier)_output_data/P_avgE_k/")
        end
        if !isdir("$(file_output_identifier)_output_data/avgE_k/")
            mkdir("$(file_output_identifier)_output_data/avgE_k/")
        end
        for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
            avg_εₖ_A_Z_TKE = DataFrame(
                A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Value = Raw_output_datafile.Avg_εₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
            )
            avg_εₖ_A = Average_over_TKE_Z(avg_εₖ_A_Z_TKE, y_A_Z_TKE)
            CSV.write(
                "$(file_output_identifier)_output_data/avgE_k/$(fissionant_nucleus_identifier)_avgE_$(k)_SCM_A_$(file_output_identifier).dat", 
                DataFrame(A = avg_εₖ_A.Argument, avg_ε = avg_εₖ_A.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            avg_εₖ_TKE = Average_over_A_Z(avg_εₖ_A_Z_TKE, y_A_Z_TKE)
            CSV.write(
                "$(file_output_identifier)_output_data/avgE_k/$(fissionant_nucleus_identifier)_avgE_$(k)_SCM_TKE_$(file_output_identifier).dat", 
                DataFrame(TKE = avg_εₖ_TKE.Argument, avg_ε = avg_εₖ_TKE.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_avg_εₖ = Probability_of_occurrence(avg_εₖ_A_Z_TKE, y_A_Z_TKE, Δavg_εₖ)
            CSV.write(
                "$(file_output_identifier)_output_data/P_avgE_k/$(fissionant_nucleus_identifier)_P_avgE_$(k)_SCM_$(file_output_identifier).dat", 
                DataFrame(avg_εₖ = probability_avg_εₖ.Argument, P = probability_avg_εₖ.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_avg_εₖ_L = Probability_of_occurrence(
                DataFrame(
                    A = avg_εₖ_A_Z_TKE.A[avg_εₖ_A_Z_TKE.A .<= A_H_min], 
                    Z = avg_εₖ_A_Z_TKE.Z[avg_εₖ_A_Z_TKE.A .<= A_H_min], 
                    TKE = avg_εₖ_A_Z_TKE.TKE[avg_εₖ_A_Z_TKE.A .<= A_H_min],
                    Value = avg_εₖ_A_Z_TKE.Value[avg_εₖ_A_Z_TKE.A .<= A_H_min]
                ), y_A_Z_TKE, Δavg_εₖ
            )
            CSV.write(
                "$(file_output_identifier)_output_data/P_avgE_k/$(fissionant_nucleus_identifier)_P_avgE_$(k)_SCM_LF_$(file_output_identifier).dat", 
                DataFrame(avg_εₖ = probability_avg_εₖ_L.Argument, P = probability_avg_εₖ_L.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_avg_εₖ_H = Probability_of_occurrence(
                DataFrame(
                    A = avg_εₖ_A_Z_TKE.A[avg_εₖ_A_Z_TKE.A .>= A_H_min], 
                    Z = avg_εₖ_A_Z_TKE.Z[avg_εₖ_A_Z_TKE.A .>= A_H_min], 
                    TKE = avg_εₖ_A_Z_TKE.TKE[avg_εₖ_A_Z_TKE.A .>= A_H_min],
                    Value = avg_εₖ_A_Z_TKE.Value[avg_εₖ_A_Z_TKE.A .>= A_H_min]
                ), y_A_Z_TKE, Δavg_εₖ
            )
            CSV.write(
                "$(file_output_identifier)_output_data/P_avgE_k/$(fissionant_nucleus_identifier)_P_avgE_$(k)_SCM_HF_$(file_output_identifier).dat", 
                DataFrame(avg_εₖ = probability_avg_εₖ_H.Argument, P = probability_avg_εₖ_H.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
        end
    end
end
if secondary_output_T
    T_A_Z_TKE = DataFrame(
        A = Raw_output_datafile.A[Raw_output_datafile.Tₖ .>= 0],
        Z = Raw_output_datafile.Z[Raw_output_datafile.Tₖ .>= 0],
        TKE = Raw_output_datafile.TKE[Raw_output_datafile.Tₖ .>= 0],
        Value = Raw_output_datafile.Tₖ[Raw_output_datafile.Tₖ .>= 0]
    )
    T_A = Average_over_TKE_Z(SeqAvg_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence,
        Value = Raw_output_datafile.Tₖ
    )), y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_T_A_$(file_output_identifier).dat", 
        DataFrame(A = T_A.Argument, T = T_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_T = Probability_of_occurrence(T_A_Z_TKE, y_A_Z_TKE, ΔT)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_T_$(file_output_identifier).dat", 
        DataFrame(T = probability_T.Argument, P = probability_T.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_T_L = Probability_of_occurrence(
        DataFrame(
            A = T_A_Z_TKE.A[T_A_Z_TKE.A .<= A_H_min], 
            Z = T_A_Z_TKE.Z[T_A_Z_TKE.A .<= A_H_min], 
            TKE = T_A_Z_TKE.TKE[T_A_Z_TKE.A .<= A_H_min],
            Value = T_A_Z_TKE.Value[T_A_Z_TKE.A .<= A_H_min]
        ), y_A_Z_TKE, ΔT
        )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_T_LF_$(file_output_identifier).dat", 
        DataFrame(T = probability_T_L.Argument, P = probability_T_L.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_T_H = Probability_of_occurrence(
        DataFrame(
            A = T_A_Z_TKE.A[T_A_Z_TKE.A .>= A_H_min], 
            Z = T_A_Z_TKE.Z[T_A_Z_TKE.A .>= A_H_min], 
            TKE = T_A_Z_TKE.TKE[T_A_Z_TKE.A .>= A_H_min],
            Value = T_A_Z_TKE.Value[T_A_Z_TKE.A .>= A_H_min]
        ), y_A_Z_TKE, ΔT
        )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_T_HF_$(file_output_identifier).dat", 
        DataFrame(T = probability_T_H.Argument, P = probability_T_H.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    if secondary_output_Eᵣ
        a = copy(Raw_output_datafile.aₖ[Raw_output_datafile.Tₖ .>= 0])
        Eᵣ_A_Z_TKE = copy(T_A_Z_TKE)
        Eᵣ_A_Z_TKE.Value .= Energy_FermiGas.(a, Eᵣ_A_Z_TKE.Value)
        probability_Eᵣ = Probability_of_occurrence(Eᵣ_A_Z_TKE, y_A_Z_TKE, ΔEᵣ)
        CSV.write(
            "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_Er_$(file_output_identifier).dat", 
            DataFrame(Eᵣ = probability_Eᵣ.Argument, P = probability_Eᵣ.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        probability_Eᵣ_L = Probability_of_occurrence(
            DataFrame(
                A = Eᵣ_A_Z_TKE.A[Eᵣ_A_Z_TKE.A .<= A_H_min], 
                Z = Eᵣ_A_Z_TKE.Z[Eᵣ_A_Z_TKE.A .<= A_H_min], 
                TKE = Eᵣ_A_Z_TKE.TKE[Eᵣ_A_Z_TKE.A .<= A_H_min],
                Value = Eᵣ_A_Z_TKE.Value[Eᵣ_A_Z_TKE.A .<= A_H_min]
            ), y_A_Z_TKE, ΔEᵣ
            )
        CSV.write(
            "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_Er_LF_$(file_output_identifier).dat", 
            DataFrame(Eᵣ = probability_Eᵣ_L.Argument, P = probability_Eᵣ_L.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        probability_Eᵣ_H = Probability_of_occurrence(
            DataFrame(
                A = Eᵣ_A_Z_TKE.A[Eᵣ_A_Z_TKE.A .>= A_H_min], 
                Z = Eᵣ_A_Z_TKE.Z[Eᵣ_A_Z_TKE.A .>= A_H_min], 
                TKE = Eᵣ_A_Z_TKE.TKE[Eᵣ_A_Z_TKE.A .>= A_H_min],
                Value = Eᵣ_A_Z_TKE.Value[Eᵣ_A_Z_TKE.A .>= A_H_min]
            ), y_A_Z_TKE, ΔEᵣ
            )
        CSV.write(
            "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_Er_HF_$(file_output_identifier).dat", 
            DataFrame(Eᵣ = probability_Eᵣ_H.Argument, P = probability_Eᵣ_H.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
    end
end
if secondary_output_avg_ε
    avg_ε_A_Z_TKE = DataFrame(
        A = Raw_output_datafile.A[Raw_output_datafile.Tₖ .>= 0],
        Z = Raw_output_datafile.Z[Raw_output_datafile.Tₖ .>= 0],
        TKE = Raw_output_datafile.TKE[Raw_output_datafile.Tₖ .>= 0],
        Value = Raw_output_datafile.Avg_εₖ[Raw_output_datafile.Tₖ .>= 0]
    )
    SeqAvg_ε_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence,
        Value = Raw_output_datafile.Avg_εₖ
    ))
    avg_ε_A = Average_over_TKE_Z(SeqAvg_ε_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_avgE_SCM_A_$(file_output_identifier).dat", 
        DataFrame(A = avg_ε_A.Argument, avg_ε = avg_ε_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    avg_ε_TKE = Average_over_A_Z(SeqAvg_ε_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_avgE_SCM_TKE_$(file_output_identifier).dat", 
        DataFrame(TKE = avg_ε_TKE.Argument, avg_ε = avg_ε_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_avg_ε = Probability_of_occurrence(avg_ε_A_Z_TKE, y_A_Z_TKE, Δavg_ε)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_avgE_SCM_$(file_output_identifier).dat", 
        DataFrame(avg_ε = probability_avg_ε.Argument, P = probability_avg_ε.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_avg_ε_L = Probability_of_occurrence(
        DataFrame(
            A = avg_ε_A_Z_TKE.A[avg_ε_A_Z_TKE.A .<= A_H_min], 
            Z = avg_ε_A_Z_TKE.Z[avg_ε_A_Z_TKE.A .<= A_H_min], 
            TKE = avg_ε_A_Z_TKE.TKE[avg_ε_A_Z_TKE.A .<= A_H_min],
            Value = avg_ε_A_Z_TKE.Value[avg_ε_A_Z_TKE.A .<= A_H_min]
        ), y_A_Z_TKE, Δavg_ε
        )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_avgE_SCM_LF_$(file_output_identifier).dat", 
        DataFrame(avg_ε = probability_avg_ε_L.Argument, P = probability_avg_ε_L.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_avg_ε_H = Probability_of_occurrence(
        DataFrame(
            A = avg_ε_A_Z_TKE.A[avg_ε_A_Z_TKE.A .>= A_H_min], 
            Z = avg_ε_A_Z_TKE.Z[avg_ε_A_Z_TKE.A .>= A_H_min], 
            TKE = avg_ε_A_Z_TKE.TKE[avg_ε_A_Z_TKE.A .>= A_H_min],
            Value = avg_ε_A_Z_TKE.Value[avg_ε_A_Z_TKE.A .>= A_H_min]
        ), y_A_Z_TKE, Δavg_ε
        )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_avgE_SCM_HF_$(file_output_identifier).dat", 
        DataFrame(avg_ε = probability_avg_ε_H.Argument, P = probability_avg_ε_H.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
end
if secondary_output_TXE_Q 
    txe_A_Z_TKE, Q_A_Z, Q_AH = Vectorized_TXE_Q(A₀, Z₀, fission_type, E_incident, y_A_Z_TKE, A_H_range, dmass_excess)
    txe_AH = Average_over_TKE_Z(txe_A_Z_TKE, y_A_Z_TKE)
    txe_TKE = Average_over_A_Z(txe_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_Q_A_Z_$(file_output_identifier).dat", 
        DataFrame(A_H = Q_A_Z.A, Z = Q_A_Z.Z, Q = Q_A_Z.Value, σQ = Q_A_Z.σ), 
        writeheader=true, newline="\r\n", delim=' '
    )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_Q_AH_$(file_output_identifier).dat", 
        DataFrame(A = Q_AH.Argument, Q = Q_AH.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_TXE_AH_$(file_output_identifier).dat", 
        DataFrame(A = txe_AH.Argument, TXE = txe_AH.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_TXE_TKE_$(file_output_identifier).dat", 
        DataFrame(TKE = txe_TKE.Argument, TXE = txe_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
end
if secondary_output_E_excitation
    E_excitation_A = Average_over_TKE_Z(E_excitation, y_A_Z_TKE)
    E_excitation_TKE = Average_over_A_Z(E_excitation, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_E_excit_A_$(file_output_identifier).dat", 
        DataFrame(A = E_excitation_A.Argument, E = E_excitation_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_E_excit_TKE_$(file_output_identifier).dat", 
        DataFrame(TKE = E_excitation_TKE.Argument, E = E_excitation_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_E_excitation = Probability_of_occurrence(E_excitation, y_A_Z_TKE, 1.0)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_P_E_excit_$(file_output_identifier).dat", 
        DataFrame(E = probability_E_excitation.Argument, P = probability_E_excitation.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    RT_A_Z_TKE, T_init_A_Z_TKE = Ratio_of_Temperatures(A₀, Z₀, A_H_range, fragmdomain, E_excitation, density_parameter_type, density_parameter_data)   
    RT_AH = Average_over_TKE_Z(RT_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_RT_AH_$(file_output_identifier).dat", 
        DataFrame(A = RT_AH.Argument, RT = RT_AH.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
end
#Write average quantities to file
open("$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_Average_quantities_$(file_output_identifier).dat", "w") do file
    if secondary_output_Yield
        avg_A_L = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .<= A_H_min])
        avg_A_H = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .>= A_H_min])
        if !isnan(avg_A_H[2])
            write(file, "<A>_L = $(avg_A_L[1]) ± $(avg_A_L[2])\n")
            write(file, "<A>_H = $(avg_A_H[1]) ± $(avg_A_H[2])\n")
        else
            write(file, "<A>_L = $(avg_A_L[1])\n")
            write(file, "<A>_H = $(avg_A_H[1])\n")
        end
        avg_TKE = Average_yield_argument(y_TKE, y_TKE.Argument)
        if !isnan(avg_TKE[2])
            write(file, "<TKE> = $(avg_TKE[1]) ± $(avg_TKE[2])\n")
        else
            write(file, "<TKE> = $(avg_TKE[1])\n")
        end
        δₑₒ = (sum(y_Z.Value[iseven.(y_Z.Argument)]) - sum(y_Z.Value[isodd.(y_Z.Argument)]))/sum(y_Z.Value)
        σδₑₒ = (1/sum(y_Z.Value)) * sqrt((1 + δₑₒ)^2 * sum(y_Z.σ .^2) + 2*δₑₒ*(sum(y_Z.σ[isodd.(y_Z.Argument)] .^2) - sum(y_Z.σ[iseven.(y_Z.Argument)].^2)))
        if !isnan(σδₑₒ)
            write(file, "δₑₒ = $(δₑₒ *100) ± $(σδₑₒ *100) %\n\n")
        else
            write(file, "δₑₒ = $(δₑₒ *100) %\n\n")
        end
    end
    if secondary_output_nu
        avg_ν_L = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_ν_H = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_ν = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_range)
        avg_ν_Pair = Average_value(ν_Pair_A_Z_TKE, y_A_Z_TKE, A_H_range)
        write(file, "<ν>_L = $avg_ν_L\n<ν>_H = $avg_ν_H\n<ν> = $avg_ν\n<ν>_pair = $avg_ν_Pair\n\n")
        if secondary_output_Ap
            avg_Ap_L = Average_yield_argument(yp_Ap, yp_Ap.Argument[yp_Ap.Argument .<= A_H_min])
            avg_Ap_H = Average_yield_argument(yp_Ap, yp_Ap.Argument[yp_Ap.Argument .>= A_H_min])
            if !isnan(avg_Ap_H[2])
                write(file, "<Ap>_L = $(avg_Ap_L[1]) ± $(avg_Ap_L[2])\n")
                write(file, "<Ap>_H = $(avg_Ap_H[1]) ± $(avg_Ap_H[2])\n\n")
            else
                write(file, "<Ap>_L = $(avg_Ap_L[1])\n")
                write(file, "<Ap>_H = $(avg_Ap_H[1])\n\n")
            end
            avg_KEp_L = Average_value(kep_A_Z_TKE, y_A_Z_TKE, A_L_range)
            avg_KEp_H = Average_value(kep_A_Z_TKE, y_A_Z_TKE, A_H_range)
            write(file, "<KEp>_L = $avg_KEp_L\n")
            write(file, "<KEp>_H = $avg_KEp_H\n\n")
        end
        if secondary_output_Tₖ
            for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
                Tₖ_A_Z_TKE = DataFrame(
                    A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Value = Raw_output_datafile.Tₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                )
                avg_Tₖ_L = Average_value(Tₖ_A_Z_TKE, y_A_Z_TKE, A_L_range)
                avg_Tₖ_H = Average_value(Tₖ_A_Z_TKE, y_A_Z_TKE, A_H_range)
                avg_Tₖ = Average_value(Tₖ_A_Z_TKE, y_A_Z_TKE, A_range)
                write(file, "<T_$(k)>_L = $avg_Tₖ_L\n<T_$(k)>_H = $avg_Tₖ_H\n<T_$(k)> = $avg_Tₖ\n\n")
            end
        end    
        if secondary_output_avg_εₖ
            for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
                avg_εₖ_A_Z_TKE = DataFrame(
                    A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Value = Raw_output_datafile.Avg_εₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                )
                avg_εₖ_L = Average_value(avg_εₖ_A_Z_TKE, y_A_Z_TKE, A_L_range)
                avg_εₖ_H = Average_value(avg_εₖ_A_Z_TKE, y_A_Z_TKE, A_H_range)
                avg_εₖ = Average_value(avg_εₖ_A_Z_TKE, y_A_Z_TKE, A_range)
                write(file, "avg_<ε_$(k)>_L = $avg_εₖ_L\navg_<ε_$(k)>_H = $avg_εₖ_H\navg_<ε_$(k)> = $avg_εₖ\n\n")
            end
        end
        if secondary_output_Eᵣ
            for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
                Tₖ = Raw_output_datafile.Tₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                aₖ = Raw_output_datafile.aₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                Eᵣ_A_Z_TKE = DataFrame(
                    A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Value = Energy_FermiGas.(aₖ, Tₖ)
                )
                avg_Eᵣ_L = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_L_range)
                avg_Eᵣ_H = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_H_range)
                avg_Eᵣ = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_range)
                write(file, "<Eᵣ_$(k)>_L = $avg_Eᵣ_L\n<Eᵣ_$(k)>_H = $avg_Eᵣ_H\n<Eᵣ_$(k)> = $avg_Eᵣ\n\n")
            end
        end
    end
    if secondary_output_T
        avg_T_L = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_T_H = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_T = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_range)
        write(file, "<T>_L = $avg_T_L\n<T>_H = $avg_T_H\n<T> = $avg_T\n\n")
    end
    if secondary_output_Eᵣ
        avg_Eᵣ_L = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_Eᵣ_H = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_Eᵣ = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_range)
        write(file, "<Eᵣ>_L = $avg_Eᵣ_L\n<Eᵣ>_H = $avg_Eᵣ_H\n<Eᵣ> = $avg_Eᵣ\n\n")
    end  
    if secondary_output_avg_ε
        avg_ε_L = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_ε_H = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_ε = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_range)
        write(file, "<ε>_L = $avg_ε_L\n<ε>_H = $avg_ε_H\n<ε> = $avg_ε\n\n")
    end  
    if secondary_output_TXE_Q 
        avg_Q = Average_value(Q_AH, y_A, A_H_range)
        avg_TXE = Average_value(txe_A_Z_TKE, y_A_Z_TKE, A_H_range)
        write(file, "<Q> = $avg_Q\n<TXE>= $avg_TXE\n\n")
    end
    if secondary_output_E_excitation
        avg_E_exi = Average_value(E_excitation, y_A_Z_TKE, A_range)
        avg_E_exi_L = Average_value(E_excitation, y_A_Z_TKE, A_L_range)
        avg_E_exi_H = Average_value(E_excitation, y_A_Z_TKE, A_H_range)
        avg_RT = Average_value(RT_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_T_init = Average_value(T_init_A_Z_TKE, y_A_Z_TKE, A_range)
        avg_T_init_L = Average_value(T_init_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_T_init_H = Average_value(T_init_A_Z_TKE, y_A_Z_TKE, A_H_range)
        write(file, "<E*>_L = $avg_E_exi_L\n<E*>_H = $avg_E_exi_H\n<E*> = $avg_E_exi\n\n")
        write(file, "<T_init>_L = $avg_T_init_L\n<T_init>_H = $avg_T_init_H\n<T_init> = $avg_T_init\n\n")
        write(file, "<RT> = $avg_RT\n\n")
    end
end