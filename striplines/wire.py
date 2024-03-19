import numpy as np

μ0 = np.pi * 4e-7

class Wire:

    def field(self, I, x, z):
         r2 = x**2 + z**2

         Bx = - μ0*I*z / (2*np.pi*r2)
         Bz = μ0*I*x / (2*np.pi*r2)

    def fourier_field(self, I, kx, z):

        '''
        use that the FT of (1/(x^2 + z^2)) is pi/z * exp(-|k|z) for z > 0   <- can see this by doing IFT of exp(-|k|z)
        for z < 0 it is - pi/z * exp(-|k|z)
        for z = 0 it is undefined I believe - should use zero in this case because Bx is zero everywhere.

        use that the FT of (x/(x^2 + z^2)) is i d/dk of the FT of (1/(x^2 + z^2)) 
        which is, for z > 0, - i * sgn(k) * pi * exp(-|k|z)
        for z < 0: it is - i * sgn(k) * pi * exp(|k|z)
        for z = 0 we will ignore issues and just treat the thing as continuous. 

        but what about where k=0? In this case we can just note this is an odd function so is zero

        '''

        Bx = np.sign(z) * - .5 * μ0 * I * np.exp(-np.abs(kx*z))
        Bz = np.sign(kx) * - .5j * μ0 * I *np.exp(-np.abs(kx*z))
        
        return Bx, Bz
    
    def integrated_fourier_field(self, I, kx, z1, z2):
        #assume that both z1 and z2 are either above or below the stripline
        pass


    
class ThickWire:

    def init(self, h):
        self.h = h

    def field(self, I, x, z):
        pass 
        #not implemented

    def fourier_field(self, I, kx, z):

        #we suppose that z is either entirely above or entirely below the stripline.

        I = I/self.h

        Bx_0 = - .5 * μ0 * I * self.h
        Bz_0 = 0

        if z >= self.h:
            #except for k = 0 
            fac = np.exp(-np.abs(kx)*z) * (np.exp(np.abs(kx)*self.h) - 1) / np.abs(kx)
            Bx =  - .5 * μ0 * I * fac
            Bz = -np.sign(kx) * .5j * μ0 * I * fac

        if z <= 0:
            fac = - np.exp(np.abs(kx)*z) * (np.exp(-np.abs(kx)*self.h) - 1) / np.abs(kx)
            Bx = .5 * μ0 * I * fac
            Bz = -np.sign(kx) * .5j * μ0 * I * fac


        Bx = np.where(kx == 0, Bx_0, Bx)
        Bz = np.where(kx ==0, Bz_0, Bz)
        
        return Bx, Bz
    
