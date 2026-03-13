import uproot
import matplotlib.pyplot as plt
import numpy as np
import sys
import pickle
import os, ROOT

# Usage: python plot_mass_pt_2D.py input.root
with open("list_ntuples_biased.txt","r") as fin:
  filenames = fin.readlines()

for filename in filenames:
  filename = filename.strip()

for filename in filenames:
  print(filename)

#filename = "AToEleEle_ntuple.root"

#allmasses = []
#allpts = []

mVSptHist = ROOT.TH2F("mVSpT", "Mass vs pT", 34, 0.01, 1.2, 38, 25, 160)


for n,filename in enumerate(filenames):
# Open ROOT file
  print(f"File {n+1} of {len(filenames)}")
  with uproot.open(filename) as f:
    # Print available keys to check structure (optional)
    # print(f.keys())
    
    # Access the tree (assuming it's named "Events" or similar)
    # You can adjust this if your file has a different tree name
      tree = f["fevt"]["RHTree"]
    

    # Read the two branches
      masses = tree["A_mass"].array(library="np")
      pts   = tree["A_pT"].array(library="np")
    
      for i in range(len(masses)):
        mass = masses[i][0]
        pt = pts[i][0]
      #  print(f"Mass: {mass}, pT: {pt}")
        #allmasses.append(mass)
        #allpts.append(pt)
        mVSptHist.Fill(mass,pt)
      

outputFile = ROOT.TFile("hist.root","RECREATE")
mVSptHist.Write()
outputFile.Close()

xlabel = "m^{a} (GeV)"
ylabel = "p_{T}^{a} (GeV)"

