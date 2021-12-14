import numpy          as     np

from   smith_kcorr    import GAMA_KCorrection
from   scipy.optimize import brentq, minimize


def rest_gmr(kcorr_rfunc, kcorr_gfunc, z, gmr):
     z = np.atleast_1d(z)
     gmr = np.atleast_1d(gmr)
     
     def pred_gmr(x):
          x = np.atleast_1d(x)
          
          # Here, x is the rest_color to be solved for in the native                                                                                    
          # reference, i.e. z=0.1 for AJS.                                                                                                              
          return  x + kcorr_gfunc(z, x) - kcorr_rfunc(z, x)
          
     def diff(x):
          result = gmr - pred_gmr(x)
          return result[0]

     def absdiff(x):
        result= np.abs(diff(x))
        return result[0]
        
     try:
        # rest color limits.  
        result = brentq(diff, -1., 2.)

     except ValueError as VE:
        # Brent method fails, requires sign change across boundaries.
        result = minimize(absdiff, 0.62)

        if result.success:
            result = result.x[0]

        else:             
            print(result.message)

            raise RuntimeError()

     return  result

def smith_rest_gmr(zs, gmrs):
   kcorr_r = GAMA_KCorrection(band='R')
   kcorr_g = GAMA_KCorrection(band='G')

   result = []
   
   for z in zs:
        for gmr in gmrs:
             interim = rest_gmr(kcorr_r.k, kcorr_g.k, z, gmr)
             
             result.append(interim)

             print(interim)

             exit(0)
             
   return np.array(result)

