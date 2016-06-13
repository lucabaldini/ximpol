import math

poldegree_file = file('Flux_output.txt','r')
outputfile = file('cyg_x1_pol_flux_wedge_corona_model_more_data.txt','w')
print "#Energy Pol degree"
#print "#Energy Flux"
#outputfile.writelines("#Energy Pol angle\n")
outputfile.writelines("#Energy Flux\n")
for line in poldegree_file:#.readlines():
    l = line.split()
    pol_degree = l[0].split('=')[1].split(',')[0]
    #log_energy =   float(l[1].split('=')[1])
    energy =   float(l[1].split('=')[1])
    #energy = 10**log_energy
    print '%s %s'%(energy, pol_degree)
    outputfile.writelines( '%s %s\n'%(energy, pol_degree))
outputfile.close()
