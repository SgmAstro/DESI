==================
DESI Change Log
==================

0.15.0 (2021-07-14)
-------------------

Note: Major changes to output formats; requires desispec >= 0.45.0

* Split FIBERMAP into FIBERMAP (coadded) and EXP_FIBERMAP (per-exposure)
  (PR `#196`_).
* Add additional ZWARN bit masking for known bad input data (PR `#196`_).
* Rename zbest -> redrock output, update rrdesi option names (PR `#198`_).

.. _`#196`: https://github.com/desihub/redrock/pull/196
.. _`#198`: https://github.com/desihub/redrock/pull/198
