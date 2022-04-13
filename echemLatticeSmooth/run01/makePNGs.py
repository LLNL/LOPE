from paraview.simple import *

view = GetActiveView()
reader = GetActiveSource()
tsteps = reader.TimestepValues

i=0
for j in tsteps:
#times = range(0,100,5) + range(100,10050,50)
#for j in times:
    #view.ViewTime = tsteps[j]
    view.ViewTime = j
    Render()
    m = str(i)
    nf = 'Image_%s.png' % m.zfill(5)
    WriteImage(nf)
    i=i+1;
