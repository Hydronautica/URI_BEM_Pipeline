.. _examples:

========
Examples
========

This chapter walks through the complete pipeline using the included
``SimplePontoon/`` test case — a twin-pontoon catamaran floating
platform at 1 : 50 laboratory scale.


Geometry
--------

* Two semicylindrical pontoons, radius :math:`R = 0.025` m, length
  :math:`L = 0.75` m
* Transverse spacing: pontoon centres at :math:`y = \pm 0.125` m
* Half-submerged (cylinder axis at waterline :math:`z = 0`)
* Mass ≈ 1.47 kg, centre of gravity at :math:`z_G = -0.005` m
* Water depth :math:`d = 1.0` m


Step 1 — Generate Mesh and Nemoh Inputs
----------------------------------------

.. code-block:: bash

   cd SimplePontoon
   python generate_case.py

This creates 1 148 nodes and 1 152 panels (including waterplane lids)
along with all Nemoh input files (``Nemoh.cal``, ``mesh/``,
``Mechanics/``).


Step 2 — Run Nemoh BEM Solver
------------------------------

.. code-block:: bash

   ../bin_windows/preProcessor
   ../bin_windows/Solver
   ../bin_windows/postProcessor

Expected run time: 2–5 minutes on a modern desktop for 120 frequencies.
Results are written to ``results/``.


Step 3 — Remove Irregular Frequencies
--------------------------------------

.. code-block:: bash

   python remove_irreg_freq.py

Before/after comparison plots are generated automatically.  Inspect the
heave added-mass curve to verify that spikes have been smoothed.

.. figure:: figures/irreg_freq_heave_detail.png
   :align: center
   :width: 75%

   Heave added mass before and after irregular-frequency removal.


Step 4 — Time-Domain Simulation
---------------------------------

Open ``run_nemoh_timeseries.m`` in MATLAB and run the script.  With the
default configuration (irregular waves, :math:`H_s=0.014` m,
:math:`T_p=3.2` s, 3-line taut mooring, 360 s duration) the solver
completes in under 30 seconds.

Key outputs:

* **Displacement, velocity, force time histories** — Figures 6–8
* **Probability density functions** — Figure 9 (transient-clipped)
* **Exceedance probability curves** — Figures 10–11, with annotated
  1 % and 5 % exceedance levels


Step 5 — Parametric Sea-State Sweep
-------------------------------------

Open ``sweep_sea_states.m`` in MATLAB and run the script.  The default
5 × 5 grid (:math:`H_s` from 0.005–0.050 m, :math:`T_p` from
0.6–3.0 s) takes approximately 10–20 minutes.

Key outputs:

* **Contour maps** — 1 % exceedance force for each active DOF
* **Summary tables** — printed to the console

Contour maps are exported to ``pdf_outputs/`` and show two quantities
per DOF:

* Net hydrodynamic force (what the mooring resists)
* True net (unbalanced) force on the body


Step 6 — Mooring Sensitivity (Optional)
-----------------------------------------

To explore different mooring configurations, modify the ``lines`` struct
in either script and re-run.  For example, switching from taut wire to
catenary:

.. code-block:: matlab

   % Replace taut-wire block with:
   for i = 1:3
       ang = (i-1)*120 * pi/180;
       lines(i).type     = 'catenary';
       lines(i).anchor   = [4*cos(ang); 4*sin(ang); -1.0];
       lines(i).fairlead = [0.15*cos(ang); 0.15*sin(ang); 0];
       lines(i).L0 = 2.5;   lines(i).w = 0.2;   lines(i).EA = 5000;
   end
   [K_mooring, F0] = compute_mooring_stiffness(lines, zeros(6,1));

Compare the resulting stiffness matrix and time-domain response to
assess the trade-off between catenary compliance and taut-wire stiffness.


Expected Directory After Full Run
-----------------------------------

.. code-block:: text

   SimplePontoon/
   ├── Nemoh.cal
   ├── input_solver.txt
   ├── generate_case.py
   ├── remove_irreg_freq.py
   ├── nemoh_to_wamit.py
   ├── plot_results.m
   ├── run_nemoh_timeseries.m
   ├── sweep_sea_states.m
   ├── compute_mooring_stiffness.m
   ├── mesh/
   │   └── Pontoon.dat
   ├── Mechanics/
   │   ├── Inertia.dat
   │   └── Kh.dat
   ├── results/
   │   ├── RadiationCoefficients.tec
   │   ├── ExcitationForce.tec
   │   ├── IRF.tec
   │   └── ...
   ├── Wamit_format/
   │   ├── Buoy.1
   │   └── Buoy.3
   └── pdf_outputs/
       ├── displacement_timeseries.pdf
       ├── force_timeseries.pdf
       ├── force_pdf.pdf
       ├── exceedance_disp.pdf
       ├── exceedance_force.pdf
       ├── sweep_hydro_*.pdf
       └── sweep_net_*.pdf
