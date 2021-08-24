from inputs import *
from band_structure import *
from real_structure import *

from hopfield_full import *
from hopfield_spherical import *

mode='full' #'full'

inputs=inputs()
#inputs.if_so()
#inputs.run_calculations()

if mode=='spherical': inputs.radwf_file='RADWF/RADWF.radwf'
band_str=band_structure()
band_str.read_ef_and_dos()
band_str.read_almblm()
if mode != 'spherical':
 band_str.read_ene()
 band_str.which_bands_cross_ef()
 band_str.calc_kweights()

#print(len(band_str.ENE[0]))
real_str=real_structure(band_str.E_f) 
real_str.read_volume()
real_str.read_potential()
real_str.read_radwf()
real_str.div_pot_and_radwf_by_r()
#exit()



if mode=='spherical': 
 band_str.read_ldos() ###ONLY FOR SPHERICAL!!!
 hop_sher=hopfield_spherical(real_str,band_str)
 hop_sher.calc_hopfield_spherical()
else:
 hop=hopfield(real_str,band_str)
 hop.calc_r_and_k_integrals()
 hop.parallel_Hopfield()
 hop.print_result()
#hop.single_Hopfield([0,2])
