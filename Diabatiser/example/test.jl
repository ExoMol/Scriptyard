include("../Diabatiser.jl")
#using ..Diabatiser # Alternative import to avoid use of `Diabatiser.` prefix
using CSV
using DataFrames
using LinearAlgebra
using PyPlot
using Tables

# Initial setup
df = CSV.read("adiabatic_pecs.csv", DataFrame)
geoms = Array(df[:, 1]) # Internuclear geometries as a 1D Array (vector)
adiabats = Array{Float64}(undef, size(df)[1], 2, 2) # 3 dimensional array of diagonal matrices representing adiabatic potential energy operator for each geometry
for i=1:size(df)[1]
    adiabats[i, :, :] = Diagonal(Array(df[i, 2:3]))
end

# Best fit parameters for NAC
opt = Diabatiser.fit_diabat(adiabats, geoms, Diabatiser.lorentzian)
println(opt)

# Diabatise from best fit params
diabats = Array{Float64}(undef, size(adiabats))
for i=1:length(geoms)
    diabats[i, :, :] = Diabatiser.diabatise(adiabats[i, :, :], geoms[i], Diabatiser.lorentzian, opt)
end

# Write results to file
table = Matrix{Float64}(undef, (length(geoms), 3))
table[:, 1] = geoms
table[:, 2] = diabats[:, 1, 1]
table[:, 3] = diabats[:, 2, 2]
CSV.write("diabatic_pecs.csv",  Tables.table(table, header=["R", "Energy 1.6", "Energy 2.6"]))

# Plot comparison of diabats and adiabats
plot(geoms, adiabats[:, 1, 1], color=:red, label="Adiabatic")
plot(geoms, adiabats[:, 2, 2], color=:red)
plot(geoms, diabats[:, 1, 1], color=:green, label="Diabatic")
plot(geoms, diabats[:, 2, 2], color=:green)
xlabel("Bond Length / Å"), ylabel("Energy / Ha"), legend(), tight_layout()
xlim(1.3, 5.0), savefig("diabatised.png")
xlim(1.4, 1.8), ylim(-149.85, -149.750), savefig("diabatised_zoom.png")
close()