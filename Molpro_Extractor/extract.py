#--- Extraction setup
states = ['1.1', '2.1', '1.7', '2.7']
quadrupoles = [
    '<1.7|QMZZ|1.7>', '<2.7|QMZZ|2.7>', '<2.7|QMZZ|1.7>', 
    '<1.7|QMYZ|2.1>', '<2.7|QMYZ|2.1>', 
    '<1.7|QMYZ|1.1>', '<2.7|QMYZ|1.1>', 
    '<1.1|QMZZ|1.1>', '<2.1|QMZZ|2.1>', '<2.1|QMZZ|1.1>']
spinorbits = ['<1.7|LSZ|2.1>']
files = [
    "b1sigp-d1pig-qtm-lx-b7 from 1.41A to 0.98A 2 states.out",
    "b1sigp-d1pig-qtm-lx-f7 from 1.41A to 3.0A 2 states.out"
]
out = "b1sigp-d1pig-qtm__2states"

#---
import numpy as np
import re
coupls = [*quadrupoles, *spinorbits]
keys = {}
vals = {}
labels = {}

#Create regex key for each MRCI state
for state in states:
    labels[state] = f"""Energy {state}"""
    keys[state] = f"""!MRCI STATE\s*{re.escape(state)}\s*Energy\s*(-?\d{{3}}\.\d{{12}})i?"""
    vals[state] = []

#Create regex key for each coupling specified
for coupl in coupls:
    labels[coupl] = coupl
    keys[coupl] = f"""!MRCI (?:expec|trans)\s*{re.escape(coupl)}\s*(-?\d{{0,3}}\.\d{{12}})i?"""
    vals[coupl] = []

#Create regex key for geometries and loop number
keys['geoms'] = """SETTING R\(\d{0,3}\)\s*=\s*(\d{1,2}.\d{8})"""
vals['geoms'] = []
keys['loop'] = """DO\sI\s*=\s*(\d{1,3}).\d{8}"""
vals['loop'] = []

for file in files:
    with open(file, "r") as f:
        while True:
            line = f.readline()
            if line == '':
                break
            for key in keys:
                #Append first regex match group to list if key is found
                match = re.search(keys[key], line)
                if match:
                    vals[key].append(match[1])
    #Check if any geometries failed
    if len(vals['geoms']) != len(vals['loop']):
        for key in keys:
            vals[key] = vals[key][:len(vals['loop'])-1]

#Spin-orbits have two values for each geometry (mean-field and Briet-Pauli)
for spinorbit in spinorbits:
    vals[spinorbit] = [[vals[spinorbit][i], vals[spinorbit][i+1]] for i in range(0, len(vals[spinorbit]), 2)]

#Write results
with open(out+".csv", "w") as f:
    #print header
    print("R", end="", file=f)
    for state in states:
        print(",", labels[state], end="", file=f)
    for coupl in quadrupoles:
        print(",", coupl, end="", file=f)
    for coupl in spinorbits:
        for val in ["MF", "BP"]:
            print(",", coupl+val, end="", file=f)
    print("", file=f)

    #print curves  
    for (g, geom) in enumerate(vals['geoms']):
        print(geom, end="", file=f)
        for state in states:
            print(",", vals[state][g], end="", file=f)
        for coupl in quadrupoles:
            print(",", vals[coupl][g], end="", file=f)
        for coupl in spinorbits:
            for val in vals[coupl][g]:
                print(",", val, end="", file=f)
        print("", file=f)
