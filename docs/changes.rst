==================
DESI Change Log
==================

5.0.2 (2022-May-20)
-------------------
* Correct for overcounting of fillfactor contribution to vmax (PR `#196`_).
* Correct for bit bug for FILLFACTOR > 0.8, with update_bit (PR `#196`_).
* Change gama gold lower z limit to 0.02 due to typo in TMR (PR `#196`_).
* Require fillfactor > 0.8 for volume of reference lf (PR `#196`_).
* Add bitmask.update_bit (PR `#196`_).
* Add __name__ clauses to prevent pool hangs (PR `#196`_).
* Add ddp_zlimits to specify external redshift limits (PR `#196`_).
* Allow GALL and RALL for multi-field catalogs (to findfile) (PR `#196`_).
* Collate multiple oversampled realizations for fillfactor calculation (in bound dist.)  (PR `#196`_)
* Start params file for common hardcoding, e.g. fillfactor threshold, sphere radius (PR `#196`_).
* External redshift limits specified by ddp_zlimits (PR `#196`_).
* Add plot idx to delta8_limits for TMR comparison (PR `#196`_).
* Notebook improvements (PR `#196`_).
* Switch to oversample 2 as default (more robust to memory & hanging issues.)  (PR `#196`_).
* Calculate exact fillfactors for galaxies, don't rely on random matching. rFILLFACTOR is matched (PR `#196`_).
* Improved initialisation of Brent method to catch color, zmin and zmax failures (PR `#196`_).
* Default jack knife as 4 jks per field (PR `#196`_).
* Mid, mean, median Ms for LF calc.  (PR `#196`_).
* Default to TMR-like LF binning (PR `#196`_).
* Add TMR d8-schechter model for plot comparison (PR `#196`_).
* Improve summary stats for more useful table comparison (PR `#196`_).
* Increase the split/complement buffer in bound_dist to 2. Mpc/h (PR `#196`_).
* Do not update zmin for fillfactor cut for VMAX evaluation in vmaxer (PR `#196`_).
* Rewrite vol. avg. fillfactor calc.  Now evaluated in vmaxer (PR `#196`_).
* <FILLFACTOR> needs to be restricted to density tiers for d8 LF (PR `#196`_).
* <FILLFACTOR> needs to be restricted to density tiers for d8 LF - field dependent correction (PR `#196`_).
* Multiple oversampling realisations and collation (PR `#196`_).
* Color-dependent stepwise (PR `#196`_).
* Slurm jobnames contain field and realisation (PR `#196`_).
* Params file for sphere radius and oversample realisation number (PR `#196`_).
* Add 'spawn' context to prevent pool hangs. TBD if effective (PR `#197`_).
* Add zmax optimisation refinements, brent with only a guess at bracketing interval and Nelder-Mead (PR `#198`_). 

.. _`#196`: https://github.com/SgmAstro/DESI/pull/196
.. _`#197`: https://github.com/SgmAstro/DESI/pull/197
.. _`#198`: https://github.com/SgmAstro/DESI/pull/198

5.0.1 (2022-April-27)
-------------------
* Restict to fillfactor > 0.8 for volfracs.
  (PR `#174`_).
* More careful header updates in gen_ddp_n8 (PR `#174`_).
* Multiple realisations (16) of DESI randoms (PR `#174`_).
* Extension to all DESI rosettes rather than GAMA (PR `#174`_).
* Change default rosette radii to allow low completeness regions (PR `#174`_).
* Remove survey specifics in favor of propagated header info (PR `#174`_).
* Limit DDP N8 counting to each single field (PR `#174`_).
* More control of propagated header info (PR `#174`_).
* Tweaked Brent initialisation to not have zmax fail on ~100 (bright) galaxies (PR `#174`_).
* Possibility of weights in multi-field luminosity function (PR `#174`_).
  
.. _`#174`: https://github.com/SgmAstro/DESI/pull/174

5.0.0 (2022-April-25)
-------------------

Note: Major changes 

* Working GitHub Actions setup and README badge
  (PR `#155`_).
* Fix lack of logs - redirect stdout in python scripts (PR `#155`_).
* Working customisation of queue in pipeline run, e.g. cosma/cordelia (PR `#155`_).
* Send pipelog scripts to the right place (PR `#155`_).
* Working configurations to write python script args to one place. replayable soon? (PR `#155`_).
* Protect against divison warnings for exactly zero fillfactor (PR `#155`_).
* Add protection against negative z for cosmo functions to remove zero div. errors (PR `#155`_).
* Add change log, this file (PR `#155`_).
* Fix fillfactor normalisation bug, where should be properly normalised against oversampled randoms (PR `#155`_).
* Comment pool.join which interferes with pytest given thread spawning (PR `#155`_).
* Findfile can find logs (PR `#155`_).
* Dryrun uses a 2x2 common patch of the field, rather than a random sample (PR `#155`_).
* Protect against zero sized array means in vol. avg. fillfactor calculation (PR `#155`_).
* Jackknife limits polish (PR `#155`_).
* Suppress merge conflict warnings (PR `#155`_).
* Fix bug in submit.py for logs (PR `#155`_).
* Fix bug in volfracs calc.: ddp1_rand = rand[rand['DDPZLIMS'][:,0] == 1] (PR `#155`_).
* Fix bug where fillfactor_vmax was incorrectly wrapped by vmaxer (PR `#155`_).
* Make selfcount_volfracs a default (PR `#155`_).
* Fix failure to pass $SURVEYARG to ddp_limits (PR `#155`_). 
* Single origin for survey field definition (PR `#155`_).
* Safe_reset: remove run files, but ignore immutable (PR `#155`_). 

.. _`#155`: https://github.com/SgmAstro/DESI/pull/155
