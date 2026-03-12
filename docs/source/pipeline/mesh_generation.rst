.. _pipeline-mesh:

===============
Mesh Generation
===============

The first stage of the pipeline creates the panel mesh that Nemoh
requires for the BEM computation.  The script ``generate_case.py``
parametrically builds a twin-pontoon (catamaran) geometry and writes all
Nemoh input files.


Script: ``generate_case.py``
----------------------------

.. rubric:: Location

``SimplePontoon/generate_case.py``


Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Description
   * - ``R``
     - 0.025 m
     - Pontoon radius
   * - ``L``
     - 0.75 m
     - Pontoon length (x-direction)
   * - ``y_c1``, ``y_c2``
     - ±0.125 m
     - Transverse centre of each pontoon
   * - ``z_c``
     - 0.0 m
     - Vertical centre (at waterline)
   * - ``rho``
     - 1000 kg/m³
     - Water density
   * - ``g``
     - 9.81 m/s²
     - Gravitational acceleration
   * - ``z_G``
     - −0.005 m
     - Centre-of-gravity height
   * - ``depth``
     - 1.0 m
     - Water depth (0 = infinite)
   * - ``N_circ``
     - 16
     - Circumferential divisions (semicircle)
   * - ``N_x``
     - 10
     - Axial divisions along pontoon
   * - ``N_r``
     - 3
     - Radial divisions on end caps
   * - ``N_lid``
     - 4
     - Transverse divisions on waterplane lid
   * - ``use_lid``
     - True
     - Add lid panels at :math:`z=0` for irregular-frequency suppression
   * - ``n_freq``
     - 120
     - Number of analysis frequencies
   * - ``omega_min / omega_max``
     - 0.2 / 40.0 rad/s
     - Frequency range


Geometry
^^^^^^^^

Each pontoon is a half-cylinder (semicircular cross-section below the
waterline) with flat end caps.  The two pontoons are mirror images about
the :math:`y=0` plane.  When ``use_lid = True`` an additional set of
flat panels is placed at the waterline (:math:`z=0`) across the interior
of each pontoon to suppress irregular frequencies in the BEM solution
(see :ref:`theory-excrad`).


Output Files
^^^^^^^^^^^^

The script creates the following in the case directory:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Content
   * - ``Nemoh.cal``
     - Master input file defining bodies, DOFs, frequency range
   * - ``mesh/Pontoon.dat``
     - Panel mesh in Nemoh DAT format (nodes + quads)
   * - ``Mechanics/Inertia.dat``
     - 6 × 6 mass/inertia matrix
   * - ``Mechanics/Kh.dat``
     - 6 × 6 hydrostatic stiffness matrix
   * - ``input_solver.txt``
     - Green's function solver settings


Running
^^^^^^^

.. code-block:: bash

   cd SimplePontoon
   python generate_case.py

The script is self-contained and requires only the Python standard
library (no NumPy).  After execution, inspect the mesh statistics
printed to the console (number of nodes, panels, displaced volume).

.. tip::

   Increasing ``N_circ`` and ``N_x`` improves mesh convergence but
   increases Nemoh run time quadratically.  A resolution study is
   recommended before production runs.
