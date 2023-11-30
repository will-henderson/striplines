import numpy as np

μ0 = np.pi * 4e-7

class StripLine:
    def __init__(self, width, centre=0, lift=0, height=0):
        self.width = width
        self.centre = centre
        self.lift = lift
        self.height = height
        
    def profile(self, x):
        '''
        Returns a function that returns 1 if x lies within the values contained by the strip, zero otherwise. 
        '''
        return TH(self.centre, self.width, x)
    
    def make_profile(self, points=1000, sideroom=1):
        '''
        Returns an array that represents the profile. 1 for values in the profile, 0 otherwise. Useful for plotting.
        
        Parameters:
            points (int): number of points in the array. 
            sideroom (number): the returned array contains sideroom zeros either side of the ones.
        '''
        x = np.linspace(centre - (sideroom+.5)*self.width, centre + (sideroom+.5)*self.width, points)
        y = self.profile(x)
        
        return x, y
    
    def field(self, I, x, z):
        x_ = x - self.centre
        z_ = z - self.lift
        w = self.width
        h = self.height
        
        if h == 0:
            Bx = - μ0 * I / (2*np.pi*w) * (np.arctan((w/2 - x_)/z_) + np.arctan((w/2 + x_)/z_))
            Bz = μ0 * I / (4*np.pi*w) * np.log(1 + 2*w*x_/((w/2 - x_)**2 + z_**2))
            
        else:
            
            Bx = - μ0 * I / (2*np.pi*w*h) * \
                ( h * (np.arctan((w/2 - x_)/(z_ - h)) + np.arctan((w/2 + x_)/(z_ - h))) \
                + z * (np.arctan((z_ - h)/(w/2 - x_)) + np.arctan((z_ - h)/(w/2 + x_))) \
                - z * (np.arctan(z_/ (w/2 - x_)) + np.arctan(z_/ (w/2 + x_))) \
                + w/4 * np.log(((w/2 + x_)**2 + z_**2 ) * ((w/2 - x_)**2 + z_**2) / (((w/2 + x_)**2 + (h - z_)**2 ) * ((w/2 - x_)**2 + (h - z_)**2))) \
                + x_/2 * np.log(((w/2 + x_)**2 + z_**2 ) * ((w/2 - x_)**2 + (h - z_)**2) / (((w/2 + x_)**2 + (h - z_)**2 ) * ((w/2 - x_)**2 + z_**2))) \
                )
            
            Bz = μ0 * I / (2*np.pi*w*h) * \
                ( h/2 * np.log(1 + 2*w*x_/((h-z_)**2 + (w/2 - x_)**2)) \
                + z_/2 * np.log((1 + 2*w*x_/(z_**2 + (w/2 - x_)**2)) / (1 + 2*w*x_/((h-z_)**2 + (w/2 - x_)**2))) \
                + (w/2+x_) * (np.arctan((h-z_)/(w/2+x_)) + np.arctan((z_)/(w/2+x_))) \
                - (w/2-x_) * (np.arctan((h-z_)/(w/2-x_)) + np.arctan((z_)/(w/2-x_))) \
                )
            
        return Bx, Bz
        
    def fourier_field(self, I, kx, z):
    
        w = self.width
        h = self.height
        z_ = -self.lift + z
        if h == 0:
            if z_ == 0:
                wirex = -μ0*I / (2*w)
            else:
                wirex = -μ0*I / (2*w) * np.exp(-np.abs(kx * z))  
            wirez = -1j * wirex * np.where(kx > 0, 1, -1)
        else: 
            wirex = -μ0*I / (2*w*h) * np.where(kx == 0, self.height, (np.exp(-np.abs(kx * (h-z_))) - np.exp(-np.abs(kx*z_))) / np.abs(kx) * np.where(z_ > 0, 1, -1))
            wirez = -1j * wirex * np.where(kx > 0, 1, -1)
            
        # then just do convolution with a top hat function
        # just product in fourier space.
        #convolved with delta function if off centre
        
        thft = w * np.sinc(w * kx/2) * np.exp(1j*kx*self.centre)
        
        Bx = wirex * thft
        Bz = wirez * thft
        
        return Bx, Bz
                



def TH(center, width, x):
    return np.where(np.abs(x - center) < width / 2, 1, 0)
