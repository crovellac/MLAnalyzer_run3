import uproot
import matplotlib.pyplot as plt
import numpy as np
import sys

import os, ROOT
import cmsstyle
import pickle

ROOT.gROOT.SetBatch(True)

file1 = ROOT.TFile("hist.root")
mVSptHist = file1.Get("mVSpT")


xlabel = "m^{a} (GeV)"
ylabel = "p_{T}^{a} (GeV)"

cmsstyle.ResetAdditionalInfo()
cmsstyle.SetExtraText("Simulation Preliminary")
cmsstyle.SetEnergy("13.6")  
cmsstyle.SetLumi("")
  
#c = cmsstyle.cmsCanvas("Testing",xmin,xmax,ymin,ymax,xlabel,"Events",extraSpace=0.06,iPos=0,yTitOffset=1.7)
c = cmsstyle.cmsCanvas("Testing",-0.13,1.34,20,180,xlabel,ylabel,square=False,iPos=0,extraSpace=0.06,yTitOffset=1)  
 
cmsstyle.cmsDraw(mVSptHist, "COLZ TEXT")

cmsstyle.SaveCanvas(c, "mVSpT.svg")

# Printing out bins for unbiasing

min_pt = 40.0
max_pt = 160.0
delta_pt = 5.0
min_mass = 0.01
max_mass = 1.2000000000000002
delta_mass = 0.07

pt_range = int((max_pt-min_pt)/delta_pt)
mass_range = int((max_mass-min_mass)/delta_mass)

pt_bin_str = "  std::vector <int> pT_bins   = {"
for p in range(1,pt_range+1):
  pt = min_pt + p*delta_pt
  pt_bin_str += str(int(pt))
  if p < pt_range:
    pt_bin_str += ", "
pt_bin_str += "};"
print(pt_bin_str)

mass_bin_str = "  std::vector <double> m_bins = {"
for m in range(1,mass_range+1):
  mass = min_mass + m*delta_mass
  mass_bin_str += str(mass)
  if m < mass_range:
    mass_bin_str += ", "
mass_bin_str += "};"
print(mass_bin_str)
print("")
print("  std::vector <int> occ = {")

comment_str = "//  "
for p in range(1,pt_range+1):
  pt = min_pt + p*delta_pt
  comment_str += str(int(pt))
  if p < pt_range:
    comment_str += ",     "
print(comment_str)


for m in range(1,mass_range+1):
  mass = min_mass + m*delta_mass
  row_str = ""
  for p in range(1,pt_range+1):
    pt = min_pt + p*delta_pt
    #print(f"Mass: {mass}\tpT: {pt}")
    thisbin = mVSptHist.FindBin(mass-delta_mass/2, pt-delta_pt/2)
    count = mVSptHist.GetBinContent(thisbin)
    if not((p == pt_range) and (m == mass_range)):
      row_str += str(int(count)) + ", "
    else:
      row_str += str(int(count)) + " "
  row_str += "// -> " + str(mass)
  print(row_str)
print("  };")
    



