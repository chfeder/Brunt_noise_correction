
'''
1. Many layers of complexity.
Just project. cg.sum(axis=2)
Include optical effects
    Resolution of CCD
    PSF of the source
    Response of the CCD, Optical Depth effects

2. Opacity
Transmission = e^( -a Column )
3. Image = Source * Transmission
Add beam width
Add CCD resolution
Add CCD response
'''
def Project():
    # Proj = CoveringGrid.sum( axis = 1 )
    # Proj = Convolve with Gausian (need the PSF framing camera)
