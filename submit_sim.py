import os

ind_init = 200
Necc = 200

for index in range(ind_init+1, ind_init+Necc+1):
    f = open('jobscript_sim', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<ind>>', '%i'%index)
    h = open('jobscript', 'w')
    h.write(g)
    h.close()

    os.system('qsub jobscript')
    os.system('rm jobscript')
    #os.system('cat jobscript')
