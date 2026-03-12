.. _theory-hydrostatics:

============
Hydrostatics
============

Hydrostatic restoring forces and the mass matrix are essential inputs to
the Cummins equation (:eq:`eq-cummins`).  This chapter describes how
they are assembled for the pipeline.


Hydrostatic Stiffness Matrix
----------------------------

For a freely floating body the linearised hydrostatic restoring matrix
:math:`\mathbf{K}_h` is a :math:`6\times6` symmetric matrix.  The
non-zero entries for a body with vertical symmetry plane are
:cite:`Faltinsen1990`:

.. math::
   :label: eq-kh

   \begin{aligned}
   K_{33} &= \rho\,g\,A_{\text{wp}}, \\[4pt]
   K_{34} = K_{43}
     &= -\rho\,g\!\iint_{A_{\text{wp}}} y\,dA, \\[4pt]
   K_{35} = K_{53}
     &= \rho\,g\!\iint_{A_{\text{wp}}} x\,dA, \\[4pt]
   K_{44} &= \rho\,g\!\iint_{A_{\text{wp}}} y^{2}\,dA
             + \rho\,g\,\nabla\,\overline{KB}
             - m\,g\,\overline{KG}, \\[4pt]
   K_{55} &= \rho\,g\!\iint_{A_{\text{wp}}} x^{2}\,dA
             + \rho\,g\,\nabla\,\overline{KB}
             - m\,g\,\overline{KG},
   \end{aligned}

where

* :math:`A_{\text{wp}}` is the waterplane area,
* :math:`\nabla` is the displaced volume,
* :math:`\overline{KB}` is the vertical distance from keel to centre of
  buoyancy,
* :math:`\overline{KG}` is the vertical distance from keel to centre of
  gravity.

Surge, sway, and yaw (:math:`K_{11}`, :math:`K_{22}`, :math:`K_{66}`)
have zero hydrostatic stiffness; restoring in those modes comes from
mooring lines (see :ref:`theory-mooring`).


Metacentric Height
------------------

The stability criterion for small rotations is captured by the
metacentric height:

.. math::
   :label: eq-gm

   \overline{GM} = \overline{KB} + \overline{BM} - \overline{KG},

where the metacentric radius is

.. math::

   \overline{BM} = \frac{I_{\text{wp}}}{\nabla},

with :math:`I_{\text{wp}}` the second moment of the waterplane area
about the relevant axis.  Positive :math:`\overline{GM}` indicates
static stability.

For the roll and pitch stiffness entries, the relationship can be
written compactly as

.. math::

   K_{44} = \rho\,g\,\nabla\,\overline{GM}_T,
   \qquad
   K_{55} = \rho\,g\,\nabla\,\overline{GM}_L,

where subscripts :math:`T` and :math:`L` denote transverse and
longitudinal metacentric heights.


Nemoh Output
------------

Nemoh computes :math:`\mathbf{K}_h` from the discretised panel mesh and
writes it to ``Mechanics/Kh.dat`` as a :math:`6\times6` matrix in
scientific notation.  The pipeline reads this file directly.


Mass Matrix
-----------

The :math:`6\times6` rigid-body mass matrix :math:`\mathbf{M}` is

.. math::
   :label: eq-mass

   \mathbf{M} =
   \begin{bmatrix}
   m & 0 & 0 & 0 & m\,z_G & -m\,y_G \\
   0 & m & 0 & -m\,z_G & 0 & m\,x_G \\
   0 & 0 & m & m\,y_G & -m\,x_G & 0 \\
   0 & -m\,z_G & m\,y_G & I_{xx} & -I_{xy} & -I_{xz} \\
   m\,z_G & 0 & -m\,x_G & -I_{xy} & I_{yy} & -I_{yz} \\
   -m\,y_G & m\,x_G & 0 & -I_{xz} & -I_{yz} & I_{zz}
   \end{bmatrix},

where :math:`(x_G,y_G,z_G)` is the centre of gravity relative to the
body-fixed origin and :math:`I_{ij}` are the moments and products of
inertia.

Nemoh writes this matrix to ``Mechanics/Inertia.dat``.


Combined Restoring Matrix
-------------------------

In the Cummins equation the total linear restoring matrix is

.. math::

   \mathbf{C} = \mathbf{K}_h + \mathbf{K}_{\text{moor}},

where :math:`\mathbf{K}_{\text{moor}}` is the linearised mooring
stiffness (see :ref:`theory-mooring`).  Both matrices are
:math:`6\times6` and symmetric.

.. seealso::

   :ref:`pipeline-timedomain` for how the mass and stiffness matrices
   are loaded from the Nemoh output files.
