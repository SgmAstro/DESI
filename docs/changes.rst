==================
DESI Change Log
==================

5.0.0 (2022-April-25)
-------------------

Note: Major changes 

* Working GitHub Actions setup and README badge
  (PR `#155`_).
* fix lack of logs - redirect stdout in python scripts (PR `#155`_).
* working customisation of queue in pipeline run, e.g. cosma/cordelia (PR `#155`_).
* send pipelog scripts to the right place (PR `#155`_).
* working configurations to write python script args to one place. replayable soon? (PR `#155`_).
* protect against divison warnings for exactly zero fillfactor.
* add protection against negative z for cosmo functions to remove zero div. errors (PR `#155`_).
* add change log, this file (PR `#155`_).
* fix fillfactor normalisation bug, where should be properly normalised against oversampled randoms (PR `#155`_).
* comment pool.join which interferes with pytest given thread spawning (PR `#155`_).
* findfile can find logs (PR `#155`_).
* dryrun uses a 2x2 common patch of the field, rather than a random sample (PR `#155`_).
* protect against zero sized array means in vol. avg. fillfactor calculation (PR `#155`_).
* jackknife limits polish (PR `#155`_).
* suppress merge conflict warnings (PR `#155`_).
* fix bug in submit.py for logs (PR `#155`_).
* fix bug in volfracs calc.: ddp1_rand = rand[rand['DDPZLIMS'][:,0] == 1] (PR `#155`_).
* 

.. _`#155`: https://github.com/desihub/redrock/pull/155
