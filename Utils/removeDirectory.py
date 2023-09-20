from ROOT import *
import numpy as np
import math

dr = "ttbar_semi"
#hh = "H1_mass"
ff = TFile("output.root","update")
ff.cd()
gDirectory.rmdir(dr)
#ff.Delete(hh)

print("Removed %s directory in output.root..."%(dr))
