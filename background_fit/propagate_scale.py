import ROOT,PlotUtils
import sys

def updatefile(f1,f2):
    keylist = f1.GetListOfKeys()
    f1.cd()
    for key in keylist:
        print(key)
        h2 = f2.Get(key.GetName())
        h1 = f1.Get(key.GetName())
        # h1.Multiply(h1,h2)
        # f1.cd()
        h1.Write("",ROOT.TObject.kOverwrite)
    f1.Write()

if __name__=="__main__":
    f1,f2 = sys.argv[1:3]
    print(f1,f2)
    updatefile(ROOT.TFile.Open(f1,"UPDATE"),ROOT.TFile.Open(f2))
