
'''
    Makes a covering grid.
    covering_grid(level, left_edge, dimensions, fields=None, ds=None, num_ghost_zones=0, use_pbar=True, field_parameters=None)
'''

def extract_section_1(ds,level,edge):
    Cg = ds.covering_grid(level=level, left_edge=edge, dimensions=ds.domain_dimensions)
    return Cg
