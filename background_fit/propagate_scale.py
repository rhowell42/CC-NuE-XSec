import ROOT,PlotUtils
import os

def updatehist(h1,h2):
    h1.Multiply(h2)
    

def updatefile(f1,f2):
    keylist = f1.GetListOfKeys()
    for key in keylist:
        h1=f1.Get(key.GetName())
        h2=f2.Get(key.GetName())
        h1.Multiply(h2)
        f1.cd()
        h1.Write("",ROOT.TObject.kOverwrite)

if __name__=="__main__":
    f1,f2 = sys.argv[1:3]
    print(f1,f2)
    updatefile(ROOT.TFile.Open(f1,"UPDATE"),ROOT.TFile.Open(f2))
