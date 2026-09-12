#Convert Light & Heavy stacked values provided on AH domain to single-file values on the whole A domain

cd(@__DIR__)
cd("raw_data/")
using CSV, DataFrames

A₀, Z₀ = NaN, NaN

rawdatafile_name = "name.dat"
rawdatafile_header = ["AH", "ZH", "valLF", "valHF"]
rawdatafile_firstline = 2
rawdatafile_delim = ' '
rawdatafile = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, 
                        skipto = rawdatafile_firstline, header = rawdatafile_header)

aux_Data = zeros(nrow(rawdatafile), 3)

for i in eachindex(aux_Data[:, 1])
    A_H = rawdatafile[!, 1][i]
    Z_H = rawdatafile[!, 2][i]
    Value_LF = rawdatafile[!, 3][i]
    Value_HF = rawdatafile[!, 4][i]
    A_L = A₀ - A_H
    Z_L = Z₀ - Z_H
    aux_Data[i, 1] = A_H
    aux_Data[i, 2] = Z_H
    aux_Data[i, 3] = Value_HF
    if A_L != A_H 
        aux_Data = vcat(aux_Data, [A_L Z_L Value_LF])
    end
end

newData = DataFrame(A = aux_Data[:, 1], Z = aux_Data[:, 2], Value = aux_Data[:, 3])
cd(@__DIR__); cd("output/")
CSV.write(rawdatafile_name, newData, delim=' ')