
'''
main.py: currently pseudo code
'''
import yt
from extract_s1.py import *


# TO OBTAIN SIGMA_RHO 
Section = extract_section_1(args)

Section1 = MotionBlur( Section )

Projection = Project( Section1 )

PDF = MakePDF( Projection )

SigmaRho = MakeSigmaRho( PDF )


# TO OBTAIN SIGMA_V 
Section_v = extract_section_2( argument )

ProbablyOtherPhysicsTBD

SigmaV = MakeSigmaV( Section_v)


# TO OBTAIN b 
B = ProduceB( SigmaRho, SigmaV)

