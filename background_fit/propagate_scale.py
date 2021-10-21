import ROOT,PlotUtils
import sys

def updatefile(f1,f2):
    keylist = ["Signal","DIS","Excess","NCCoh"]
    f1.cd()
    for key in keylist:
        print(key)
        h2 = f2.Get(key)
        h1 = f1.Get(key)
        h1.AddMissingErrorBandsAndFillWithCV(h2)
        h2.AddMissingErrorBandsAndFillWithCV(h1)
        h1.Multiply(h1,h2)
        # f1.cd()
        h1.Write("",ROOT.TObject.kOverwrite)

if __name__=="__main__":
    f1,f2 = sys.argv[1:3]
    print(f1,f2)
    updatefile(ROOT.TFile.Open(f1,"UPDATE"),ROOT.TFile.Open(f2))
