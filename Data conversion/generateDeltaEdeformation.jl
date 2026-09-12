#Generate ΔE_deformation file for given β₀ & β_scission parametrizations

cd(@__DIR__); cd("raw_data/")
using CSV, DataFrames

dβ₀ = CSV.read("B2MOLLER.ANA", DataFrame; delim=' ', ignorerepeated=true, header=["Z", "A", "Value"], skipto = 2)

function Linear_function(x, a, b)
    return a*x + b
end
function Beta_sciss(Z) #NaN stopgaps beyond interval limits
    if Z < 28
        return NaN
    elseif Z >= 28 && Z <= 41
        a = 0.58/13
        b = -28*a
        return Linear_function(Z, a, b)
    elseif Z >= 42 && Z <= 44
        return 0.58
    elseif Z >= 45 && Z <= 50
        a = -0.58/6
        b = -50*a
        return Linear_function(Z, a, b)
    elseif Z >= 51 && Z <= 70
        a = 0.6/15
        b = -50*a
        return Linear_function(Z, a, b)
    else
        return NaN
    end
end
function E_LDM(β, A, Z)
    η = (A - 2*Z)/A
    χ = 1 - 1.7826 * η^2
    α² = (5 * β^2)/(4*π)
    X_νs = -χ*(15.4941*A - 17.9439*A^(2/3) * (1 + 0.4*α²))
    X_e = Z^2 *(0.7053 * (1 - 0.2*α²)/(A^(1/3)) - 1.1529/A)
    return X_νs + X_e
end
function ΔE_def(dβ₀)
    ΔE = DataFrame(A = Int[], Z = Int[], Value = Float64[])
    for A in sort(unique(dβ₀.A))
        for Z in sort(unique(dβ₀.Z[dβ₀.A .== A]))
            β_sciss = Beta_sciss(Z)
            β_0 = dβ₀.Value[(dβ₀.A .== A) .& (dβ₀.Z .== Z)][1]
            if !isnan(β_sciss)
                E_LDM_0 = E_LDM(0, A, Z)
                E_def_sciz = E_LDM(β_sciss, A, Z) - E_LDM_0
                E_def_at = E_LDM(β_0, A, Z) - E_LDM_0
                δE_def = abs(E_def_sciz - E_def_at)
                push!(ΔE.A, A)
                push!(ΔE.Z, Z)
                push!(ΔE.Value, δE_def)
            end
        end
    end
    return ΔE
end

ΔE_deformation = ΔE_def(dβ₀)
newData = DataFrame(A = ΔE_deformation.A, Z = ΔE_deformation.Z, deltaE_def = ΔE_deformation.Value)
cd(@__DIR__); cd("output/")
newdatafile_name = "name.dat"
CSV.write(newdatafile_name, newData, delim=' ')
