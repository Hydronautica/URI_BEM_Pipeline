.. _theory-spectra:

=====================
Wave Spectra
=====================

The pipeline supports both regular (monochromatic) and irregular
(spectral) wave inputs.  This chapter describes the spectral models and
the :math:`k`-space discretisation strategy used for irregular-wave
synthesis.


JONSWAP / Pierson–Moskowitz Spectrum
------------------------------------

The JONSWAP spectrum :cite:`Hasselmann1973` is

.. math::
   :label: eq-jonswap

   S(\omega) = \frac{\alpha_{\text{PM}}\,g^{2}}{\omega^{5}}
   \exp\!\Bigl[-\frac{5}{4}\Bigl(\frac{\omega_p}{\omega}\Bigr)^{\!4}\Bigr]
   \;\gamma^{\,\exp\!\bigl[-\frac{(\omega-\omega_p)^{2}}
   {2\,\sigma_s^{2}\,\omega_p^{2}}\bigr]},

where

* :math:`\omega_p = 2\pi / T_p` is the spectral peak frequency,
* :math:`\gamma` is the peak-enhancement factor (JONSWAP default ≈ 3.3;
  :math:`\gamma = 1` recovers the Pierson–Moskowitz spectrum),
* :math:`\sigma_s` equals 0.07 for :math:`\omega\le\omega_p` and 0.09
  otherwise,
* :math:`\alpha_{\text{PM}}` is calibrated so that
  :math:`H_s = 4\sqrt{m_0}` with :math:`m_0 = \int S\,d\omega`.

The pipeline calibrates :math:`\alpha_{\text{PM}}` by computing the
zeroth moment numerically and rescaling.


Dispersion Relation
-------------------

In finite water depth :math:`d`, wave frequency and wavenumber are
linked by :cite:`Dean1991`

.. math::
   :label: eq-dispersion

   \omega^{2} = g\,k\,\tanh(k\,d).

Given :math:`\omega`, the wavenumber :math:`k` is found iteratively with
a Newton–Raphson solver:

.. math::

   k^{(n+1)} = k^{(n)}
   - \frac{g\,k^{(n)}\tanh(k^{(n)}d) - \omega^{2}}
          {g\bigl[\tanh(k^{(n)}d) + k^{(n)}d\,\mathrm{sech}^{2}(k^{(n)}d)\bigr]},

converging in 4–6 iterations to machine precision with the deep-water
initial guess :math:`k^{(0)}=\omega^{2}/g`.


Wavenumber-Space Discretisation
-------------------------------

A naïve uniform discretisation in :math:`\omega`-space produces
**commensurate** frequencies whose superposition repeats with period
:math:`T_{\text{rep}}=2\pi/\Delta\omega`.  To avoid artificial
periodicity without requiring thousands of components, the pipeline
discretises uniformly in **wavenumber** :math:`k`-space:

1. Determine the wavenumber bounds from the user-specified frequency
   window :math:`[\omega_{\min},\,\omega_{\max}]`:

   .. math::

      k_{\min} = k(\omega_{\min}), \qquad
      k_{\max} = k(\omega_{\max}).

2. Divide :math:`[k_{\min},\,k_{\max}]` into :math:`N` uniform bins of
   width :math:`\Delta k = (k_{\max}-k_{\min})/N`:

   .. math::

      k_n = k_{\min} + \Bigl(n - \tfrac{1}{2}\Bigr)\Delta k,
      \qquad n = 1,\dots,N.

3. Map each :math:`k_n` back to frequency via the dispersion relation
   :eq:`eq-dispersion`:

   .. math::

      \omega_n = \sqrt{g\,k_n\,\tanh(k_n\,d)}.

   Because the mapping :math:`k\mapsto\omega` is nonlinear, the
   resulting :math:`\{\omega_n\}` are **incommensurate** — their ratios
   are generally irrational, so the synthesised signal never repeats
   exactly.

4. Convert the spectral bin width from :math:`k`-space to
   :math:`\omega`-space using the group velocity:

   .. math::
      :label: eq-cg

      C_g = \frac{\omega}{2k}
      \Bigl(1 + \frac{2kd}{\sinh 2kd}\Bigr),

   giving :math:`\Delta\omega_n = C_{g,n}\,\Delta k`.

5. Compute the wave amplitude for each component:

   .. math::
      :label: eq-amp

      a_n = \sqrt{2\,S(\omega_n)\,\Delta\omega_n}.


With :math:`N=128` components this approach provides statistically
converged results over simulation durations of 30 minutes or more,
whereas an equivalent :math:`\omega`-space discretisation would require
:math:`N>1000` to avoid visible repetition artefacts.


Excitation Force Reconstruction
-------------------------------

The time-domain excitation force on DOF :math:`i` is reconstructed from
the frequency-domain transfer function :math:`X_i(\omega)` (magnitude
and phase from Nemoh) and the random-phase wave components:

.. math::
   :label: eq-fexc

   F_{\text{exc},i}(t)
   = \sum_{n=1}^{N} a_n\,|X_i(\omega_n)|
     \cos\!\bigl(\omega_n\,t + \varphi_n + \angle X_i(\omega_n)\bigr),

where :math:`\varphi_n` are independent uniform random phases on
:math:`[0,2\pi)`.

.. seealso::

   :ref:`pipeline-timedomain` for the MATLAB configuration of wave
   parameters, and :ref:`theory-bem` for how the excitation transfer
   functions are computed.
