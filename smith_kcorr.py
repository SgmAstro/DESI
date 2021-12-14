class GAMA_KCorrection(object):
    def __init__(self, k_corr_file, kind="linear"):

        """
        Colour-dependent polynomial fit to the GAMA K-correction (Fig. 13 of Smith+17), 
        used to convert between SDSS r-band Petrosian apparent magnitudes, and rest 
        frame absolute manigutues at z_ref = 0.1
        
        Args:
            k_corr_file: file of polynomial coefficients for each colour bin
            z0: reference redshift. Default value is z0=0.1
            kind: type of interpolation between colour bins,
                  e.g. "linear", "cubic". Default is "linear"
        """
        
        # read file of parameters of polynomial fit to k-correction
        # polynomial k-correction is of the form
        # A*(z-z0)^4 + B*(z-z0)^3 + C*(z-z0)^2 + D*(z-z0) + E 
        col_min, col_max, A, B, C, D, E, col_med = \
            np.loadtxt(k_corr_file, unpack=True)
    
        z0 = 0.1
        self.z0 = z0              # reference redshift
        self.nbins = len(col_min) # number of colour bins in file
        
        self.colour_min = np.min(col_med)
        self.colour_max = np.max(col_med)
        self.colour_med = col_med

        # functions for interpolating polynomial coefficients in rest-frame color.
        self.__A_interpolator = self.__initialize_parameter_interpolator(A, col_med, kind=kind)
        self.__B_interpolator = self.__initialize_parameter_interpolator(B, col_med, kind=kind)
        self.__C_interpolator = self.__initialize_parameter_interpolator(C, col_med, kind=kind)
        self.__D_interpolator = self.__initialize_parameter_interpolator(D, col_med, kind=kind)
        self.__E = E[0]

        # Linear extrapolation for z > 0.5
        self.__X_interpolator = lambda x: None
        self.__Y_interpolator = lambda x: None
        self.__X_interpolator, self.__Y_interpolator = self.__initialize_line_interpolators() 

   
    def __initialize_parameter_interpolator(self, parameter, median_colour, kind="linear"):
        # returns function for interpolating polynomial coefficients, as a function of colour
        return interp1d(median_colour, parameter, kind=kind, fill_value="extrapolate")

    
    def __initialize_line_interpolators(self):
        # linear coefficients for z>0.5
        X = np.zeros(self.nbins)
        Y = np.zeros(self.nbins)
        
        # find X, Y at each colour
        redshift = np.array([0.48,0.5])
        arr_ones = np.ones(len(redshift))
        for i in range(self.nbins):
            if i <= 3:
                Q0 = 2.12
            else:
                Q0 = 0.8
            
            k = self.k(redshift, arr_ones*self.colour_med[i]) - Q0*(redshift-self.z0)
            X[i] = (k[1]-k[0]) / (redshift[1]-redshift[0])
            Y[i] = k[0] - X[i]*redshift[0]
        
        X_interpolator = interp1d(self.colour_med, X, kind='linear', fill_value="extrapolate")
        Y_interpolator = interp1d(self.colour_med, Y, kind='linear', fill_value="extrapolate")
        
        return X_interpolator, Y_interpolator


    def __A(self, colour):
        # coefficient of the z**4 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__A_interpolator(colour_clipped)

    def __B(self, colour):
        # coefficient of the z**3 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__B_interpolator(colour_clipped)

    def __C(self, colour):
        # coefficient of the z**2 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__C_interpolator(colour_clipped)

    def __D(self, colour):
        # coefficient of the z**1 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__D_interpolator(colour_clipped)

    def __X(self, colour):
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__X_interpolator(colour_clipped)

    def __Y(self, colour):
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__Y_interpolator(colour_clipped)


    def k(self, redshift, colour):
        """
        Polynomial fit to the GAMA K-correction for z<0.5
        The K-correction is extrapolated linearly for z>0.5

        Args:
            redshift: array of redshifts
            colour:   array of ^0.1(g-r) colour
        Returns:
            array of K-corrections
        """
        K = np.zeros(len(redshift))
        idx = redshift <= 0.5
        
        K[idx] = self.__A(colour[idx])*(redshift[idx]-self.z0)**4 + \
                 self.__B(colour[idx])*(redshift[idx]-self.z0)**3 + \
                 self.__C(colour[idx])*(redshift[idx]-self.z0)**2 + \
                 self.__D(colour[idx])*(redshift[idx]-self.z0) + self.__E 

        idx = redshift > 0.5
        
        K[idx] = self.__X(colour[idx])*redshift[idx] + self.__Y(colour[idx])
        
        if (colour[0] > 0.76):
            Q0 = 0.80
        else:
            Q0 = 2.12

        return  K + Q0*(redshift-self.z0)

def test_plots():

    kcorr_r = GAMA_KCorrection("k_corr_rband_z01.mpeg")
    kcorr_g = GAMA_KCorrection("k_corr_gband_z01.mpeg")

    z = np.arange(-0.01,0.601,0.01)
    cols = 0.130634, 0.298124, 0.443336, 0.603434, 0.784644, 0.933226, 1.06731
    
    # make r-band k-correction plot
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_r.k(z, col)
        plt.plot(z, k, label=r"$^{0.1}(g-r)_\mathrm{med}=%.3f$"%c)

    plt.xlabel(r"$z$")
    plt.ylabel(r"$^{0.1}K_r(z)+E(z)$")
    plt.xlim(0,0.6)
    plt.ylim(-0.6,1)
    plt.legend(loc="upper left").draw_frame(False)
    plt.show()

    # make g-band k-correction plot
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_g.k(z, col)
        plt.plot(z, k, label=r"$^{0.1}(g-r)_\mathrm{med}=%.3f$"%c)

    plt.xlabel(r"$z$")
    plt.ylabel(r"$^{0.1}K_g(z)+E(z)$")
    plt.xlim(-0.01,0.6)
    plt.ylim(-0.4,1.4)
    plt.legend(loc="upper left").draw_frame(False)
    plt.show()
        
if __name__ == "__main__":
    test_plots()
    

gama_kcorr_r = GAMA_KCorrection("k_corr_rband_z01.mpeg")
gama_kcorr_g = GAMA_KCorrection("k_corr_gband_z01.mpeg")