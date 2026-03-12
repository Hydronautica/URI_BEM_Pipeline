.. _theory-mooring:

=======
Mooring
=======

The pipeline includes an analytical mooring-line solver that computes
line profiles, tensions, and a linearised :math:`6\times6` stiffness
matrix :math:`\mathbf{K}_{\text{moor}}` for use in the Cummins equation
:cite:`Faltinsen1990,DNVOSE301`.


Catenary Lines
--------------

For a mooring line of unstretched length :math:`L_0`, weight per unit
length :math:`w`, and maximum seabed-contact length, the classical
catenary equations relate the horizontal :math:`(D)` and vertical
:math:`(H)` span to the fairlead tensions :math:`F_H` (horizontal) and
:math:`F_V` (vertical) :cite:`Irvine1981`:

**Case 1 — Partial seabed contact** (:math:`F_V < w\,L_0`):

The suspended length is :math:`L_s = F_V / w` and the resting length is
:math:`L_b = L_0 - L_s`.

.. math::
   :label: eq-cat-partial

   \begin{aligned}
   D &= L_b + \frac{F_H}{w}\ln\!\Bigl(
        \frac{F_V + \sqrt{F_H^{2}+F_V^{2}}}{F_H}\Bigr), \\[6pt]
   H &= \frac{1}{w}\Bigl(\sqrt{F_H^{2}+F_V^{2}} - F_H\Bigr).
   \end{aligned}

**Case 2 — Fully suspended** (:math:`F_V \ge w\,L_0`):

.. math::
   :label: eq-cat-full

   \begin{aligned}
   D &= \frac{F_H}{w}\Bigl[
        \sinh^{-1}\!\Bigl(\frac{F_V}{F_H}\Bigr)
        - \sinh^{-1}\!\Bigl(\frac{F_V - w\,L_0}{F_H}\Bigr)
        \Bigr], \\[6pt]
   H &= \frac{1}{w}\Bigl[
        \sqrt{F_H^{2}+F_V^{2}}
        - \sqrt{F_H^{2}+(F_V-w\,L_0)^{2}}
        \Bigr].
   \end{aligned}

The solver uses Newton–Raphson iteration with analytic Jacobian
:math:`\partial(D,H)/\partial(F_H,F_V)` to invert these equations for
the tensions given the fairlead position.


Taut-Wire (Elastic) Lines
--------------------------

For taut synthetic-fibre or steel-wire lines the line hangs nearly
straight and the dominant restoring mechanism is axial elasticity rather
than geometry.  The tension–extension relationship is

.. math::
   :label: eq-taut

   T = EA\,\frac{L - L_0}{L_0},

where

* :math:`L = |\mathbf{r}_{\text{fair}} - \mathbf{r}_{\text{anchor}}|`
  is the current chord length,
* :math:`L_0` is the natural (unstretched) length,
* :math:`EA` is the axial stiffness.

The force vector at the fairlead is directed along the line:

.. math::

   \mathbf{F}_{\text{line}} = -T\,\hat{\mathbf{e}},
   \qquad
   \hat{\mathbf{e}} = \frac{\mathbf{r}_{\text{fair}}
   - \mathbf{r}_{\text{anchor}}}{L}.


Numerical Linearisation
-----------------------

The full :math:`6\times6` mooring stiffness matrix
:math:`\mathbf{K}_{\text{moor}}` is obtained by central-difference
numerical differentiation of the total mooring force about the
equilibrium position :cite:`Jonkman2009`:

.. math::
   :label: eq-kmoor

   K_{ij} = -\frac{F_i(\mathbf{q}_0+\delta\mathbf{e}_j)
   - F_i(\mathbf{q}_0-\delta\mathbf{e}_j)}{2\,\delta},

where :math:`F_i` is the :math:`i`-th component of the net mooring
force, :math:`\mathbf{e}_j` is the :math:`j`-th unit vector in the
generalised-displacement space, and :math:`\delta` is a small
perturbation (translations in metres, rotations in radians).

For a multi-line configuration the contributions from each line are
summed:

.. math::

   \mathbf{K}_{\text{moor}}
   = \sum_{\ell=1}^{N_{\text{lines}}}
     \mathbf{K}_{\text{moor}}^{(\ell)}.

The ``compute_mooring_stiffness.m`` function handles both catenary and
taut-wire lines and automatically selects the appropriate equations
based on the ``type`` field of each line structure.


Implementation Notes
--------------------

* The Newton–Raphson catenary solver converges in 5–15 iterations with a
  tolerance of :math:`10^{-10}`.
* Fairlead coordinates are transformed from body-fixed to global
  coordinates at each perturbation step using a small-angle rotation
  matrix.
* The perturbation size :math:`\delta` defaults to 1 mm for
  translations and :math:`10^{-4}` rad for rotations.
* A plotting option visualises the 3-D mooring layout and individual
  line profiles.

.. seealso::

   :ref:`pipeline-mooring` for practical usage of the mooring solver,
   and :ref:`theory-hydrostatics` for how
   :math:`\mathbf{K}_{\text{moor}}` is combined with
   :math:`\mathbf{K}_h`.
