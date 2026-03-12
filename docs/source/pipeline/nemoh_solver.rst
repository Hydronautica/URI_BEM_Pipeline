.. _pipeline-nemoh:

============
Nemoh Solver
============

After mesh generation the Nemoh BEM solver computes the frequency-domain
hydrodynamic coefficients.  This page describes how to run the solver
and the post-processing scripts that clean up the output.


Nemoh Three-Stage Run
---------------------

Nemoh is executed as three consecutive programs, all located in
``bin_windows/`` (or compiled from source on Linux/macOS):

1. **preProcessor** — reads ``Nemoh.cal`` and the mesh file, assembles
   the influence-coefficient matrices.
2. **Solver** — solves the BEM system at each frequency for diffraction
   and radiation potentials.
3. **postProcessor** — integrates pressures to obtain forces, computes
   impulse-response functions and infinite-frequency added mass.

.. code-block:: bash

   cd SimplePontoon
   ../bin_windows/preProcessor
   ../bin_windows/Solver
   ../bin_windows/postProcessor


Output Files
------------

All results are written to the ``results/`` subdirectory in Tecplot
``.tec`` format:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Content
   * - ``RadiationCoefficients.tec``
     - :math:`A_{ij}(\omega)` and :math:`B_{ij}(\omega)` for all DOF pairs
   * - ``ExcitationForce.tec``
     - Magnitude and phase of :math:`F_{\text{exc},i}(\omega)`
   * - ``IRF.tec``
     - Impulse-response functions :math:`K_{ij}(t)` and :math:`A_{ij,\infty}`
   * - ``DiffractionForce.tec``
     - Diffraction-only component of excitation
   * - ``FKForce.tec``
     - Froude–Krylov component of excitation

See :ref:`theory-excrad` for the data format and column definitions.


Irregular-Frequency Removal
----------------------------

Script: ``remove_irreg_freq.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Even with lid panels the BEM solution may contain residual
irregular-frequency artefacts.  The post-processing script detects and
smooths these spikes.

.. rubric:: Algorithm

1. Apply a Savitzky–Golay filter (window = 15, order = 3) to obtain a
   smooth reference curve.
2. Compute the residual between raw and filtered data.
3. Flag points where the residual exceeds ``CLIP_SIGMA`` × MAD (median
   absolute deviation).
4. Replace flagged points with the filtered values.

.. rubric:: Configuration

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``SG_WINDOW``
     - 15
     - Savitzky–Golay window length (odd integer)
   * - ``SG_POLY``
     - 3
     - Polynomial order
   * - ``CLIP_SIGMA``
     - 4
     - MAD multiplier for outlier threshold

.. rubric:: Running

.. code-block:: bash

   python remove_irreg_freq.py

The script reads from ``results/``, generates before/after comparison
plots, and writes cleaned coefficient files.


WAMIT-Format Conversion
------------------------

Script: ``nemoh_to_wamit.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For interoperability with tools that expect WAMIT-format output, the
``nemoh_to_wamit.py`` script converts Nemoh results to ``.1`` (added
mass / damping) and ``.3`` (excitation force) files, with
non-dimensionalisation by :math:`\rho`, :math:`g`, and a reference
length.

.. code-block:: bash

   python nemoh_to_wamit.py

Output is written to ``Wamit_format/``.

.. seealso::

   :ref:`pipeline-timedomain` for how the MATLAB solver reads the Nemoh
   output directly.
