from fermio import bxsf, frmsf
import numpy as np


surface = bxsf.from_file('test/FERMISURFACE.bxsf')
center = surface.center
axis = surface.axis
k_list = surface.k_list
E_list = surface.E_list
EF = surface.EF
band_num = surface.band_num
index = surface.index

# #print dimensions of the data
# print(f'axis: {axis.shape}')
# print(f'k_list: {k_list.shape}')
# print(f'E_list: {E_list.shape}')
# print(f'band_num: {len(band_num)}')
# print(f'index: {index}')
# print(f'EF: {EF}')
# print(f'center: {center}')

frmsf_file = 'test/FERMISURFACE.frmsf'
frmsf_surface = frmsf.from_file(frmsf_file)
frmsf_surface.to_bxsf('test/test_to_bxsf.bxsf')
