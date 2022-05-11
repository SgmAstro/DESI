==================
DESI Change Log
==================

5.0.1 (2022-April-27)
-------------------
* Restict to fillfactor > 0.8 for volfracs.
  (PR `#165`_).
* More careful header updates in gen_ddp_n8.
* Multiple realisations (16) of DESI randoms (PR #174)
* Extension to all DESI rosettes rather than GAMA (PR #174)
* Change default rosette radii to allow low completeness regions (PR #174)
* Remove survey specifics in favor of propagated header info (PR #174)
* Limit DDP N8 counting to each single field (PR #174)
* More control of propagated header info (PR #174)
* Tweaked Brent initialisation to not have zmax fail on ~100 (bright) galaxies (PR #174)
* Possibility of weights in multi-field luminosity function (PR #174)
  
.. _`#165`: https://github.com/desihub/redrock/pull/165

5.0.0 (2022-April-25)
-------------------

Note: Major changes 

* Working GitHub Actions setup and README badge
  (PR `#155`_).
* Fix lack of logs - redirect stdout in python scripts (PR `#155`_).
* Working customisation of queue in pipeline run, e.g. cosma/cordelia (PR `#155`_).
* Send pipelog scripts to the right place (PR `#155`_).
* Working configurations to write python script args to one place. replayable soon? (PR `#155`_).
* Protect against divison warnings for exactly zero fillfactor.
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

.. _`#155`: https://github.com/desihub/redrock/pull/155
