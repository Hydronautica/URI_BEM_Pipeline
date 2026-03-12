.. _theory-cummins:

================
Cummins Equation
================

The Cummins equation :cite:`Cummins1962` transforms the frequency-domain
hydrodynamic coefficients produced by the BEM solver into a causal,
time-domain equation of motion suitable for transient and irregular-wave
simulations.


Equation of Motion
------------------

For a floating body with :math:`N` degrees of freedom the time-domain
equation of motion is

.. math::
   :label: eq-cummins

   \bigl[\mathbf{M} + \mathbf{A}_{\infty}\bigr]\,\ddot{\mathbf{q}}(t)
   + \int_{0}^{t} \mathbf{K}(t-\tau)\,\dot{\mathbf{q}}(\tau)\,d\tau
   + \mathbf{C}\,\mathbf{q}(t)
   = \mathbf{F}_{\text{exc}}(t)
   + \mathbf{F}_{\text{ext}}(t),

where

* :math:`\mathbf{M}` is the :math:`6\times6` rigid-body mass/inertia
  matrix.
* :math:`\mathbf{A}_{\infty} = \lim_{\omega\to\infty}\mathbf{A}(\omega)`
  is the infinite-frequency added-mass matrix.
* :math:`\mathbf{K}(t)` is the matrix of radiation impulse-response
  functions (retardation kernels).
* :math:`\mathbf{C} = \mathbf{K}_h + \mathbf{K}_{\text{moor}}` is the
  combined hydrostatic and mooring restoring matrix.
* :math:`\mathbf{F}_{\text{exc}}(t)` is the wave excitation force
  reconstructed from frequency-domain transfer functions.
* :math:`\mathbf{F}_{\text{ext}}(t)` captures any additional external
  forces (wind, current, PTO, etc.).


Impulse-Response Functions
--------------------------

The retardation kernel is related to the frequency-domain radiation
damping by the cosine transform :cite:`Cummins1962`:

.. math::
   :label: eq-irf

   K_{ij}(t) = \frac{2}{\pi}\int_{0}^{\infty}
   B_{ij}(\omega)\cos(\omega t)\,d\omega,

and the infinite-frequency added mass can be recovered as

.. math::
   :label: eq-ainf

   A_{ij,\infty} = A_{ij}(\omega)
   + \frac{1}{\omega}\int_{0}^{\infty}
     K_{ij}(t)\sin(\omega t)\,dt
   \qquad \forall\;\omega>0.

In practice Nemoh directly outputs :math:`K_{ij}(t)` and
:math:`A_{ij,\infty}` in the ``IRF.tec`` results file, so the user does
not need to evaluate these integrals manually.


State-Space Form
----------------

Direct numerical evaluation of the convolution integral in
:eq:`eq-cummins` is expensive because it requires storing the full
velocity history.  Instead, each kernel :math:`K_{ij}(t)` is
approximated by a Prony series (see :ref:`theory-prony`) which converts
the convolution into a set of ordinary differential equations.

Defining the radiation state vector
:math:`\mathbf{I}_{ij}(t)\in\mathbb{R}^{N_p}` for the :math:`(i,j)`
kernel, the convolution is replaced by

.. math::
   :label: eq-ss-conv

   \int_{0}^{t} K_{ij}(t-\tau)\,\dot{q}_j(\tau)\,d\tau
   \;\approx\; \sum_{k=1}^{N_p} I_{ij,k}(t),

with each Prony state evolving as

.. math::
   :label: eq-ss-ode

   \dot{I}_{ij,k}(t) = \beta_{ij,k}\,I_{ij,k}(t)
   + \alpha_{ij,k}\,\dot{q}_j(t),
   \qquad I_{ij,k}(0) = 0.


Complete ODE System
-------------------

Collecting the body displacement :math:`\mathbf{q}`, velocity
:math:`\dot{\mathbf{q}}`, and all Prony states
:math:`\mathbf{I}` into a single state vector, the equation of motion
becomes a first-order ODE suitable for standard integrators such as
MATLAB's ``ode45``:

.. math::
   :label: eq-ode-system

   \begin{cases}
   \dot{\mathbf{q}} = \mathbf{v}, \\[4pt]
   [\mathbf{M}+\mathbf{A}_{\infty}]\,\dot{\mathbf{v}}
     = \mathbf{F}_{\text{exc}}
     + \mathbf{F}_{\text{ext}}
     - \mathbf{C}\,\mathbf{q}
     - \displaystyle\sum_{i,j,k} I_{ij,k}, \\[4pt]
   \dot{I}_{ij,k} = \beta_{ij,k}\,I_{ij,k}
     + \alpha_{ij,k}\,v_j.
   \end{cases}


DOF Selection
-------------

Not all six degrees of freedom may be of interest (e.g. a spar buoy may
be constrained to heave only).  The pipeline allows the user to specify
an ``active_dofs`` vector, and the matrices :math:`\mathbf{M}`,
:math:`\mathbf{A}_{\infty}`, :math:`\mathbf{C}`, and
:math:`\mathbf{K}(t)` are sub-indexed accordingly before integration.

.. seealso::

   :ref:`theory-prony` for the Prony fitting procedure, and
   :ref:`pipeline-timedomain` for practical usage of the time-domain
   solver.
