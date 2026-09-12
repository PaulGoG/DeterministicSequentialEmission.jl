cd(@__DIR__)
cd("input_data/")
using CSV, DataFrames, Plots

#Old output format idea
#=
    open("output_data/$(fissionant_nucleus_identifier)_main_DSE_.OUT", "w") do file
        for A in unique(Processed_raw_output.A)
            for Z in unique(Processed_raw_output.Z[Processed_raw_output.A .== A])
                P_Z_A = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
                Q = Q_value_released(A_0, Z_0, A, Z, dm)
                a = density_parameter(density_parameter_type, A, Z, density_parameter_data)
                S = Separation_energy(1, 0, A, Z, dm)[1]
                write(file, "Fission fragment: A = $A / Z = $Z / p(Z,A) = $P_Z_A / Q = $(Q[1]) / a = $a / Sₙ = $S\r\n\r\n")
                for TKE in unique(Processed_raw_output.TKE[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z)])
                    TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, E_CN[1], E_CN[2])[1]
                    E_excit = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    write(file, "TKE = $TKE / TXE = $TXE / E* = $E_excit\r\n")
                    for k in unique(Processed_raw_output.No_Sequence[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE)])
                        if k > 0
                            write(file, "   *Emission sequence $k:\r\n   ")
                            T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            write(file, "T = $T_k ")
                            a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            write(file, "/ a = $a_k ")
                            Eᵣ_k = Energy_FermiGas(a_k, T_k)
                            write(file, "/ Eʳ = $Eᵣ_k ")
                            S_k = Separation_energy(1, 0, A-k, Z, dm)[1]
                            write(file, "/ Sₙ = $S_k ")
                            avgε_k = Processed_raw_output.Avg_εₖ[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            write(file, "/ <ε> = $avgε_k\r\n")
                        else
                            write(file, "   *Fragment does not emit neutrons at TKE = $(TKE)!\r\n")
                        end
                    end
                    write(file, "\r\n")
                end
                write(file, "$horizontal_delimiter\r\n")
            end
        end
        write(file, horizontal_delimiter)
    end
=#
#####

