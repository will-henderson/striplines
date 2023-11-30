import numpy as np
from strip import StripLine

class CPW:
    def __init__(self, signal_width, ground_width, separation, lift=0, height=0):
        self.signal_width = signal_width
        self.ground_width = ground_width
        self.separation = separation
        self.centre_sep = separation + signal/2 + ground/2
        self.lift = lift
        self.height = height
        self.total_width = signal_width + 2*(ground_width + separation)
        
        self.ground_left = StripLine(ground_width, -self.centre_sep, lift, height)
        self.ground_right = StripLine(ground_width, self.centre_sep, lift, height)
        self.signal = StripLine(signal_width, 0, lift, height)
        
    def profile(self, x):
        
        return self.ground_left.profile(x) + self.ground_right.profile(x) + self.signal.profile(x)
    
    def make_profile(self, points=1000, sideroom=1):
        x = np.linspace(-(sideroom + .5) * self.total_width, (sideroom + .5) * self.total_width, points)
        y = self.profile(x)
        
        return x, y
    
    def field(self, I, x, z):
        
        Bx_gl, Bz_gl = self.ground_left.field(-.5 * I, x, z)
        Bx_gr, Bz_gr = self.ground_right.field(-.5 * I, x, z)
        Bx_s,  Bz_s  = self.signal.field(I, x, z)
        
        Bx = Bx_gl + Bx_gr + Bx_s
        Bz = Bz_gl + Bz_gr + Bz_s
        
        return Bx, Bz
    
    def fourier_field(self, I, kx, z):
        
        Bx_gl, Bz_gl = self.ground_left.fourier_field(-.5 * I, x, z)
        Bx_gr, Bz_gr = self.ground_right.fourier_field(-.5 * I, x, z)
        Bx_s,  Bz_s  = self.signal.field(I, x, z)
        
        Bx = Bx_gl + Bx_gr + Bx_s
        Bz = Bz_gl + Bz_gr + Bz_s
        
        return Bx, Bz

    
        
        