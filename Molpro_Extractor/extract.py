#--- Extraction setup
states = ['1.4', '1.6', '2.6']
quadrupoles = [
    '<1.4|QMZZ|1.4>', '<1.6|QMZZ|1.6>', '<2.6|QMZZ|2.6>', '<2.6|QMZZ|1.6>' 
]
spinorbits = ['<1.4|LSX|1.6>', '<1.4|LSX|2.6>']
files = [
    "x3sig-d1pig6b-lsx from 1.35A to 0.98A 2states.out",
    "x3sig-d1pig6f-lsx from 1.35A to 3.0A 2states.out"
]
out = "x3sig-d1pig-ls__2states"

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
            if key in spinorbits:
                vals[key] = vals[key][:2*len(vals['loop'])-1]    
            else:
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
