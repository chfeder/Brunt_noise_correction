
'''
    Makes a covering grid.
    covering_grid(level, left_edge, dimensions, fields=None, ds=None, num_ghost_zones=0, use_pbar=True, field_parameters=None)
'''
from main import *

#class extract_section_1():

def extract_section_1(ds,level,edge,size):
    Cg = ds.covering_grid(level, edge, size)
    
    return Cg
