include("Diabatiser.jl")
using DelimitedFiles

c_array = readdlm("C.txt") # load PEC data and internuclear data
r_c = c_array[:,1]				   # take column 1 as internuclear sep (Angstroms)
C3Pi_a_pec = c_array[:,2]          # columns 2 as C3Pi adiabat PEC
Cp3Pi_a_pec = c_array[:,3]         # columns 3 as C'3Pi adiabat PEC
zero = [0 for i=1:length(r_c)]     # define zeros array 

## Construct 2x2 adiabatic PEC matrices for each internuclear geometry

adiabats = Vector{Matrix}(undef, length(r_c))                    # initialize the 'array'

for i=1:length(r_c)
	adiabats[i]=[C3Pi_a_pec[i] zero[i];zero[i] Cp3Pi_a_pec[i]]
end

## calculate optimal NAC parameters

opt = Diabatiser.fit_diabat(r_c, adiabats, Diabatiser.lorentzian)

## Now calculate the diabatic curves

diabats = Diabatiser.diabatise.(r_c, adiabats, Diabatiser.lorentzian, [opt for i=1:length(adiabats)])

