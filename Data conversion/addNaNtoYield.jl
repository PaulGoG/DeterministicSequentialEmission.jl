#Add erroY column filled with NaN values to Yield distributions that have no said column

cd(@__DIR__); cd("raw_data/")
using CSV, DataFrames

rawdatafile_name = "Yield.Yield"
rawdatafile_header = ["A", "TKE", "Val"]
rawdatafile_firstline = 2
rawdatafile_delim = ' '
rawdatafile = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, 
                        skipto = rawdatafile_firstline, header = rawdatafile_header)

aux_col = [NaN for i in 1:nrow(rawdatafile)]

newYield = DataFrame(A = rawdatafile.A, TKE = rawdatafile.TKE, Value = rawdatafile.Val, erroY = aux_col)
cd(@__DIR__); cd("output/")
CSV.write(rawdatafile_name, newYield, delim=' ')