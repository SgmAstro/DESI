import os
import sys
import json
import runtime
import argparse
import pylab as pl
import numpy as np
import astropy.io.fits as fits

from   astropy.table    import Table, vstack
from   vmaxer           import vmaxer, vmaxer_rand
from   lumfn            import lumfn
from   lumfn_stepwise   import lumfn_stepwise
from   schechter        import schechter, named_schechter
from   renormalise_d8LF import renormalise_d8LF
from   delta8_limits    import d8_limits
from   config           import Configuration
from   findfile         import findfile, fetch_fields, overwrite_check, gather_cat, call_signature
from   jackknife_limits import solve_jackknife, set_jackknife

def process_cat(fpath, vmax_opath, field=None, survey='gama', rand_paths=[], extra_cols=[], bitmasks=[], fillfactor=False, conservative=False, stepwise=False, version='GAMA4'):        
    assert 'vmax' in vmax_opath

    opath = vmax_opath

    if not os.path.isfile(fpath):
        # Do not crash and burn, but proceed on gracefully. 
        print('WARNING:  Failed to find {}'.format(fpath))
        return  1

    zmax = Table.read(fpath)

    if len(zmax) == 0:
        print('Zero length catalogue, nothing to be done; Exiting.') 
        return 0
         
    found_fields = np.unique(zmax['FIELD'].data)
        
    print('Found fields: {}'.format(found_fields))
    
    minz = zmax['ZSURV'].min()
    maxz = zmax['ZSURV'].max()
    
    print('Found redshift limits: {:.3f} < z < {:.3f}'.format(minz, maxz))

    if field != None:
        assert  len(found_fields) == 1, 'ERROR: expected single-field restricted input, e.g. G9.'

    vmax  = vmaxer(zmax, minz, maxz, fillfactor=fillfactor, conservative=conservative, extra_cols=extra_cols)

    print('WARNING:  Found {:.3f}% with zmax < 0.0'.format(100. * np.mean(vmax['ZMAX'] <= 0.0)))

    vmax.meta['EXTNAME'] = 'VMAX'
    # vmax.meta['INPUT_CAT'] = fpath.replace(os.environ['GOLD_DIR'], '$GOLD_DIR')
        
    print('Writing {}.'.format(opath))

    vmax.write(opath, format='fits', overwrite=True)
    
    ##  Luminosity fn.
    opath  = opath.replace('vmax', 'lumfn')
    result = lumfn(vmax, bitmask='IN_D8LUMFN')

    result.meta['EXTNAME'] = 'LUMFN'
    # result.meta['INPUT_CAT'] = fpath.replace(os.environ['GOLD_DIR'], '$GOLD_DIR')

    result.write(opath, format='fits', overwrite=True)
    
    return  0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold luminosity function.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
    parser.add_argument('--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', default='G9')
    parser.add_argument('--survey', help='Select survey', default='gama')
    parser.add_argument('--density_split', help='Trigger density split luminosity function.', action='store_true')
    parser.add_argument('--dryrun', action='store_true', help='dryrun.')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--selfcount_volfracs', help='Apply volfrac corrections based on randoms counting themselves as ddps.', action='store_true')
    parser.add_argument('--version', help='Version', default='GAMA4')
    parser.add_argument('--conservative', help='Conservative analysis choices', action='store_true')
    
    args          = parser.parse_args()

    log           = args.log
    field         = args.field.upper()
    dryrun        = args.dryrun
    survey        = args.survey
    density_split = args.density_split
    self_count    = args.selfcount_volfracs
    version       = args.version
    conservative  = args.conservative
    
    if not density_split:
        if log:
            logfile = findfile(ftype='lumfn', dryrun=False, survey=survey, log=True)
            
            print(f'Logging to {logfile}')
                
            sys.stdout = open(logfile, 'w')

        print('Generating Gold reference LF.')

        call_signature(dryrun, sys.argv)

        # Bounded by gama gold, reference schechter limits:  
        # 0.039 < z < 0.263.
        # Note: not split by field. 

        prefix = 'randoms'
        
        # MJW/HACK:  repeated calls in this script to specify version == GAMA4? 
        fpath  = findfile(ftype='ddp',  dryrun=dryrun, survey=survey, prefix=prefix, version=version)
        opath  = findfile(ftype='vmax', dryrun=dryrun, survey=survey, prefix=prefix, version=version)

        if args.nooverwrite:
            overwrite_check(opath)

        print(f'Reading: {fpath}')
        print(f'Writing: {opath}')

        process_cat(fpath, opath, survey=survey, fillfactor=False)

        vmax            = Table.read(opath)
        rand_vmax       = vmaxer_rand(survey=survey, ftype='randoms_bd_ddp_n8', dryrun=dryrun, prefix=prefix, conservative=conservative)

        njack, jk_volfrac, limits, jks = solve_jackknife(rand_vmax)

        rand_vmax['JK']               = jks
        rand_vmax.meta['NJACK']       = njack
        rand_vmax.meta['JK_VOLFRAC']  = jk_volfrac

        vmax['JK']                    = set_jackknife(vmax['RA'], vmax['DEC'], limits=limits, debug=False)
        vmax.meta['NJACK']            = njack
        vmax.meta['JK_VOLFRAC']       = jk_volfrac

        jpath                         = findfile(ftype='jackknife', prefix=prefix, dryrun=dryrun)

        with open(jpath, 'w') as ofile:
            json.dump(limits, ofile)

        print(f'Writing: {jpath}')

        lpath                         = findfile(ftype='lumfn', dryrun=dryrun, survey=survey, prefix=prefix, version=version)

        lumfn(vmax, jackknife=np.arange(njack), opath=lpath)

        print(f'Written {lpath}')

        print('Done.')

        if log:
            sys.stdout.close()

    else:
        if log:
            # TODO NOTE: Do not support version.
            logfile = findfile(ftype='ddp_n8_d0_vmax', dryrun=False, field=field, survey=survey, log=True).replace('vmax', 'lumfn').replace('_{utier}', '')
                        
            print(f'Logging to {logfile}')
        
            sys.stdout = open(logfile, 'w')

        print('Generating Gold density-split LF.')

        call_signature(dryrun, sys.argv)

        assert  field != None

        prefix = 'randoms_ddp1'

        rpath = findfile(ftype='randoms_bd_ddp_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix, version=version)
        
        if dryrun:
            # A few galaxies have a high probability to be in highest density only. 
            utiers = np.array([8])

        else:
            utiers = np.arange(len(d8_limits))
                    
        for idx in utiers:
            ddp_idx   = idx + 1

            # Bounded by DDP1 z limits. 
            ddp_fpath = findfile(ftype='ddp_n8_d0', dryrun=dryrun, field=field, survey=survey, utier=idx, prefix=prefix, version=version)
            ddp_opath = findfile(ftype='ddp_n8_d0_vmax', dryrun=dryrun, field=field, survey=survey, utier=idx, prefix=prefix, version=version)
    
            print()
            print('Reading: {}'.format(ddp_fpath))

            try:
                failure   = process_cat(ddp_fpath, ddp_opath, field=field, rand_paths=[rpath], extra_cols=['MCOLOR_0P0', 'FIELD'], fillfactor=True, stepwise=False)

            except Exception as E:
                print('Error: Failed gen_gold_lf --density_split on d0 tier {:d} with Exception:'.format(idx))
                print(E)
                print('skipping.')
                
                continue 
        
            print('LF process cat. complete.')
                    
            result    = Table.read(ddp_opath.replace('vmax', 'lumfn'))        
            # result.pprint()

            # Single-field values.
            rand      = Table.read(findfile(ftype='randoms_bd_ddp_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix))

            fdelta    = float(rand.meta['DDP1_d{}_VOLFRAC'.format(idx)])
            fdelta_zp = float(rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)])

            result.meta['DDP1_d{}_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta)
            result.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta_zp)

            # MJW:  Load three-field randoms/meta directly. 
            rand_vmax = vmaxer_rand(survey=survey, ftype='randoms_bd_ddp_n8', dryrun=dryrun, prefix=prefix, conservative=conservative)
            
            njack, jk_volfrac, limits, jks = solve_jackknife(rand_vmax)

            rand_vmax['JK'] = jks
            rand_vmax.meta['NJACK'] = njack
            rand_vmax.meta['JK_VOLFRAC'] = jk_volfrac
            
            jpath = findfile(ftype='jackknife', prefix=prefix, dryrun=dryrun)
            
            # to do: write jk_limits to file
            
            with open(jpath, 'w') as ofile:
                json.dump(dictionary, ofile)
            
            
            
            fdelta    = float(rand_vmax.meta['DDP1_d{}_VOLFRAC'.format(idx)])
            fdelta_zp = float(rand_vmax.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)])

            d8        = float(rand_vmax.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(idx)])
            d8_zp     = float(rand_vmax.meta['DDP1_d{}_TIERMEDd8'.format(idx)])

            if (fdelta > 0.0) & (fdelta_zp > 0.0):
                result    = renormalise_d8LF(idx, result, fdelta, fdelta_zp, self_count)
            
            else:
                assert dryrun, 'ERROR:  lf renormalisation has failed.'

            result['REF_SCHECHTER']  = named_schechter(result['MEDIAN_M'], named_type='TMR')
            result['REF_SCHECHTER'] *= (1. + d8) / (1. + 0.007)

            result['REF_RATIO']      = result['PHI_IVMAX'] / result['REF_SCHECHTER']

            print('LF renormalization and ref. schechter complete.')
            
            result.pprint()

            # Reference Schechter - finer binning
            sch_Ms = np.arange(-23., -15., 1.e-3)

            sch    = named_schechter(sch_Ms, named_type='TMR')
            sch   *= (1. + d8) / (1. + 0.007)

            ##
            ref_result = Table(np.c_[sch_Ms, sch], names=['MS', 'REFSCHECHTER'])            
            ref_result.meta['DDP1_d{}_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta)
            ref_result.meta['DDP1_d{}_TIERMEDd8'.format(idx)] = '{:.6e}'.format(d8)
            ref_result.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta_zp)
            ref_result.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(idx)] = '{:.6e}'.format(d8_zp)
            
            print('Writing {}'.format(ddp_opath.replace('vmax', 'lumfn')))

            ##  
            keys           = sorted(result.meta.keys())
            
            header         = {}
            
            for key in keys:
                header[key] = str(result.meta[key])

            primary_hdu    = fits.PrimaryHDU()
            hdr            = fits.Header(header)
            result_hdu     = fits.BinTableHDU(result, name='LUMFN', header=hdr)
            ref_result_hdu = fits.BinTableHDU(ref_result, name='REFERENCE')
            hdul           = fits.HDUList([primary_hdu, result_hdu, ref_result_hdu])

            hdul.writeto(ddp_opath.replace('vmax', 'lumfn'), overwrite=True, checksum=True)

        print('Done.')

        if log:
            sys.stdout.close()

