.. _theory-prony:

============
Prony Method
============

The Prony method converts each radiation impulse-response function (IRF)
:math:`K_{ij}(t)` into a sum of complex exponentials, enabling the
convolution in the Cummins equation (:eq:`eq-cummins`) to be replaced by
a small system of ODEs :cite:`PerezFossen2009`.


Prony Series Representation
----------------------------

Each IRF is approximated as

.. math::
   :label: eq-prony

   K_{ij}(t) \approx \sum_{k=1}^{N_p}
   \alpha_{ij,k}\,e^{\,\beta_{ij,k}\,t},
   \qquad t \ge 0,

where

* :math:`\alpha_{ij,k}\in\mathbb{R}` are the **residues** (amplitudes),
* :math:`\beta_{ij,k}\in\mathbb{R}_{<0}` are the **poles** (decay
  rates), constrained to be strictly negative so that the kernel decays
  to zero as :math:`t\to\infty`,
* :math:`N_p` is the number of Prony terms per kernel (typically
  3–8 per DOF pair).


Fixed-Pole Fitting
------------------

The pipeline uses a **fixed-pole** variant of the Prony method:

1. A set of candidate poles :math:`\{\beta_k\}` is generated, uniformly
   spaced on a logarithmic scale between user-specified bounds
   (``prony_beta_min``, ``prony_beta_max``).

2. For each candidate pole set the residues :math:`\{\alpha_k\}` are
   obtained by solving the linear least-squares problem

   .. math::
      :label: eq-lsq

      \min_{\boldsymbol{\alpha}}
      \sum_{m=1}^{N_t}
      \Bigl[K_{ij}(t_m)
      - \sum_{k=1}^{N_p}\alpha_k\,e^{\,\beta_k\,t_m}\Bigr]^{2}.

3. Terms whose residues are below a threshold (relative to the maximum)
   are discarded, yielding a sparse, stable approximation.

This approach avoids the nonlinear root-finding of the classical Prony
algorithm and guarantees stability by construction.


State-Space Conversion
----------------------

Each retained Prony term generates one scalar ODE.  Defining

.. math::

   I_{ij,k}(t) = \alpha_{ij,k}\!\int_{0}^{t}
   e^{\,\beta_{ij,k}(t-\tau)}\,\dot{q}_j(\tau)\,d\tau,

differentiation gives

.. math::
   :label: eq-prony-ode

   \dot{I}_{ij,k} = \beta_{ij,k}\,I_{ij,k}
   + \alpha_{ij,k}\,\dot{q}_j,
   \qquad I_{ij,k}(0) = 0,

and the full radiation convolution is

.. math::

   \int_{0}^{t} K_{ij}(t-\tau)\,\dot{q}_j(\tau)\,d\tau
   \approx \sum_{k=1}^{N_p} I_{ij,k}(t).


Fit Quality Metric
------------------

The quality of the Prony fit is assessed with the normalised root-mean
square error (NRMSE):

.. math::
   :label: eq-nrmse

   \text{NRMSE}_{ij} = \frac{
     \sqrt{\dfrac{1}{N_t}\sum_{m=1}^{N_t}
       \bigl[K_{ij}(t_m) - \hat{K}_{ij}(t_m)\bigr]^{2}}}
     {\max_m |K_{ij}(t_m)|},

where :math:`\hat{K}_{ij}` is the Prony reconstruction.  The pipeline
prints the NRMSE for every active DOF pair and warns when any value
exceeds a user-defined tolerance (default 5 %).


Practical Recommendations
-------------------------

* **Number of terms** (``prony_nterms``): start with 4–6 per DOF pair;
  increase if NRMSE is too large.
* **Pole range**: the slowest pole should capture the longest decay
  time-scale of :math:`K(t)`, while the fastest pole should be several
  times the highest frequency of interest.
* **Threshold**: a relative threshold of :math:`10^{-3}` typically
  removes noise terms without degrading accuracy.

.. seealso::

   :cite:`Yu1995` for the state-space modelling framework and
   :ref:`pipeline-timedomain` for configuration of the Prony parameters
   in the MATLAB solver.
