#Concatenate data for ΔZ(A) & rms(A) in a single file

cd(@__DIR__); cd("raw_data/")
using CSV, DataFrames

rawdatafile_name_1 = "file1.dat"
rawdatafile_header = ["A", "Val"]
rawdatafile_firstline = 2
rawdatafile_delim = ' '
rawdatafile_1 = CSV.read(rawdatafile_name_1, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, 
                        skipto = rawdatafile_firstline, header = rawdatafile_header)

rawdatafile_name_2 = "file2.dat"
rawdatafile_2 = CSV.read(rawdatafile_name_2, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, 
                        skipto = rawdatafile_firstline, header = rawdatafile_header)

newdata =  DataFrame(A = rawdatafile_1.A, deltaZ = rawdatafile_1.Val, rms = rawdatafile_2.Val)
newdatafile_name = "DeltaZ_rms_A.nuclid"
cd(@__DIR__); cd("output/")
CSV.write(newdatafile_name, newdata, delim=' ')

