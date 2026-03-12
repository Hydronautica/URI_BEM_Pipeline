.. _introduction:

============
Introduction
============

Overview
--------

The **URI BEM Analysis Pipeline** is an integrated toolchain for the
preliminary hydrodynamic analysis of floating structures.  It couples the
open-source boundary element method (BEM) solver
`Nemoh <https://gitlab.com/lheea/Nemoh>`_ with custom MATLAB and Python
scripts to deliver a complete workflow:

1. **Mesh generation** — parametric panel-mesh creation for arbitrary
   multi-body geometries (Python).
2. **Frequency-domain BEM solve** — computation of added mass, radiation
   damping, and wave excitation coefficients via Nemoh, with automated
   irregular-frequency removal (Python).
3. **Time-domain simulation** — Cummins-equation integration with Prony
   state-space radiation approximation, irregular-wave excitation
   (JONSWAP / Pierson–Moskowitz), hydrostatic restoring, mooring forces,
   and user-defined external loads (MATLAB).
4. **Parametric sea-state sweeps** — batch evaluation of 1 % exceedance
   forces across an :math:`(H_s,\,T_p)` grid for structural design
   screening (MATLAB).
5. **Mooring stiffness analysis** — catenary and taut-wire line solvers
   with numerical linearisation to a 6 × 6 stiffness matrix (MATLAB).

.. figure:: figures/pipeline_flowchart.png
   :align: center
   :width: 85%
   :name: fig-pipeline

   High-level pipeline flowchart.


Motivation
----------

Early-stage design of floating offshore structures requires rapid
evaluation of hydrodynamic loads across many candidate geometries and sea
states.  Commercial BEM codes and their post-processing suites are
powerful but expensive and opaque.  This pipeline provides a
**fully open-source, scriptable alternative** that is:

* **Transparent** — every equation is documented and every intermediate
  result is accessible.
* **Reproducible** — deterministic seeds, version-controlled inputs, and
  automated batch runs.
* **Extensible** — the modular architecture makes it straightforward to
  swap solvers, add nonlinear effects, or couple with structural FEA.


Target Application
------------------

The pipeline was developed at the **University of Rhode Island (URI)** for
a scale-model twin-pontoon catamaran floating platform.  While the
included example (``SimplePontoon/``) reflects that geometry, the tools
are general and applicable to any floating body that can be meshed with
flat quadrilateral panels.


How to Cite
-----------

If you use this pipeline in academic work, please cite Nemoh and the
underlying theory:

* Nemoh BEM solver — :cite:`Babarit2015`
* Cummins equation — :cite:`Cummins1962`
* State-space identification — :cite:`PerezFossen2009`

A BibTeX entry for the pipeline itself will be provided once a DOI is
assigned.


License
-------

This project is released under the MIT License.  See ``LICENSE`` in the
repository root for details.
