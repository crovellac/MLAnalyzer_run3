import h5py
import numpy as np
from ROOT import TH2D, TCanvas, TFile, gROOT, TH2F
import os
#import cmsstyle

gROOT.SetBatch(True)

files = [
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_1.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_2.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_3.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_4.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_5.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_6.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_7.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_8.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_9.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_10.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_11.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_12.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_13.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_14.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_15.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_16.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_17.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_18.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_19.h5",
"AToEleEle_m0p01to1p2_pT40to160_dr0p4_ctau0_unbiased_20.h5"
]


mVSptHist = TH2F("mVSpT", "Mass vs pT", 17, 0.01, 1.2, 16, 40, 160)

total_count = 0

nevents = 17400
for n,filename in enumerate(files):
    idx = 0
    with h5py.File(filename, "r") as f:
        print(f"File {n+1} of {len(files)}")
        while idx < nevents:
            if (idx % 10 == 0):
                print(f"Image {idx} of {nevents}")
            data = list(f["all_jet"][idx:idx+1])
            pts = list(f["apt"][idx:idx+1])
            masses = list(f["am"][idx:idx+1])
            #print(f"len(masses) = {len(masses)}")
            for i in range(len(masses)):
                mVSptHist.Fill(masses[i][0], pts[i][0])
            idx = idx + 1

mVSptHist.SaveAs("hist.root")

#print(f"{total_count} total events")

#xlabel = "m^{a} (GeV)"
#ylabel = "p_{T}^{a} (GeV)"

#cmsstyle.ResetAdditionalInfo()
#cmsstyle.SetExtraText("Simulation Preliminary")
#cmsstyle.SetEnergy("13.6")  
#cmsstyle.SetLumi("")
  
#c = cmsstyle.cmsCanvas("Testing",xmin,xmax,ymin,ymax,xlabel,"Events",extraSpace=0.06,iPos=0,yTitOffset=1.7)
#c = cmsstyle.cmsCanvas("Testing",-0.01,1.21,20,180,xlabel,ylabel,square=False,iPos=0,extraSpace=0.06,yTitOffset=1)  
 
#cmsstyle.cmsDraw(mVSptHist, "COLZ TEXT")

#cmsstyle.SaveCanvas(c, "mVSpT_img.svg")
#cmsstyle.SaveCanvas(c, "mVSpT_img.root")
