.. _theory-bem:

================================
Boundary Element Method (BEM)
================================

This chapter summarises the linear potential-flow boundary-value problem
(BVP) solved by `Nemoh <https://gitlab.com/lheea/Nemoh>`_
:cite:`Babarit2015`.

Governing Equations
-------------------

Under the assumptions of an inviscid, incompressible, irrotational fluid
the velocity field derives from a potential :math:`\Phi(\mathbf{x},t)`.
For time-harmonic motions at frequency :math:`\omega` we write

.. math::
   :label: eq-decomposition

   \Phi(\mathbf{x},t)
   = \operatorname{Re}\!\bigl[\phi(\mathbf{x})\,e^{-i\omega t}\bigr],

where the complex spatial potential :math:`\phi` satisfies the Laplace
equation throughout the fluid domain :math:`\Omega`:

.. math::
   :label: eq-laplace

   \nabla^{2}\phi = 0
   \qquad \text{in } \Omega.


Boundary Conditions
-------------------

The BVP is closed by conditions on every boundary of the domain:

**Free surface** (:math:`z=0`):

.. math::
   :label: eq-fs

   -\omega^{2}\,\phi + g\,\frac{\partial\phi}{\partial z} = 0.

**Sea bed** (:math:`z=-d`):

.. math::
   :label: eq-bed

   \frac{\partial\phi}{\partial z} = 0.

**Body surface** (:math:`S_B`):

.. math::
   :label: eq-body

   \frac{\partial\phi}{\partial n} = v_n,

where :math:`v_n` is the prescribed normal velocity (zero for the
diffraction problem, :math:`-i\omega\,\eta_j\,n_j` for the radiation
problem of mode :math:`j`).

**Radiation condition** (as :math:`r\to\infty`):

.. math::
   :label: eq-radiation

   \sqrt{r}\left(\frac{\partial\phi}{\partial r}
   - ik\,\phi\right) \to 0,

with :math:`k` the fundamental wavenumber from the dispersion relation
:eq:`eq-dispersion`.


Potential Decomposition
-----------------------

The total potential is split into incident, diffracted, and radiated
components:

.. math::
   :label: eq-phitot

   \phi = \phi_I + \phi_D + \sum_{j=1}^{6} \dot{\eta}_j\,\phi_{R,j}.

Nemoh solves separately for the diffraction potential :math:`\phi_D` and
each unit-velocity radiation potential :math:`\phi_{R,j}`.


Panel Method
------------

The boundary integral equation is discretised by approximating the body
surface :math:`S_B` with :math:`N_p` flat quadrilateral panels, each
carrying a constant source strength :math:`\sigma_k`.  Applying the
body boundary condition at each panel centroid yields a dense
:math:`N_p \times N_p` linear system whose coefficient matrix involves
the free-surface Green's function :math:`G(\mathbf{x},\boldsymbol{\xi})`
:cite:`Newman1977`:

.. math::
   :label: eq-bie

   \frac{1}{2}\,\sigma(\mathbf{x})
   + \iint_{S_B} \sigma(\boldsymbol{\xi})\,
     \frac{\partial G}{\partial n_{\mathbf{x}}}
     \,dS_{\boldsymbol{\xi}}
   = v_n(\mathbf{x}).

The Green's function satisfies the free-surface, sea-bed, and radiation
conditions analytically, so only the body surface needs to be meshed.


Hydrodynamic Coefficients
-------------------------

Once the potentials are known on each panel the pressure is obtained
from the linearised Bernoulli equation and integrated over :math:`S_B`
to yield:

* **Added mass** :math:`A_{ij}(\omega)` — in-phase force per unit
  acceleration.
* **Radiation damping** :math:`B_{ij}(\omega)` — out-of-phase force per
  unit velocity.
* **Excitation force** :math:`F_{\text{exc},i}(\omega)` — complex
  amplitude of wave-induced force, combining Froude–Krylov and
  diffraction contributions.

See :ref:`theory-excrad` for detailed definitions.


Irregular Frequencies
---------------------

At certain discrete frequencies the interior Dirichlet problem has
non-trivial solutions, causing spurious spikes in :math:`A(\omega)` and
:math:`B(\omega)`.  The pipeline includes a post-processing script
(``remove_irreg_freq.py``) that detects and removes these artefacts by
interpolation; see :ref:`theory-excrad` and :ref:`pipeline-nemoh`.

.. seealso::

   :cite:`Malenica1998`, :cite:`Lee1988` for the theoretical background
   on irregular frequencies.
