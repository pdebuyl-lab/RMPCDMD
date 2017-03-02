.. _quaternions:

Quaternion velocity Verlet
==========================

Intro and notation
------------------

Rigid-body dynamics can be performed either by applied distance constraints to the bonds in
the body with RATTLE or by using a "quaternion integrator" :cite:`allen_tildesley_1987`. In
the quaternion approach, the rigid body is completely determined by its center-of-mass
position and velocity, by its quaternion orientation and by its angular velocity.

Here, we follow the velocity Verlet form of the algorithm proposed by Rozmanov and Kusalik
:cite:`rozmanov_rotational_2010`.

For clarity, the notation for quaternions and the rigid body's coordinates is given explicitly.

- :math:`r` is a vector denoting the position of a point or particle in Euclidean space.
- :math:`r_i` is as  :math:`r`, but for the  :math:`i^\mathrm{th}` particle of the rigid body.
- :math:`r_{i,j}` is the  :math:`k^\mathrm{th}` component of :math:`r_i`.
- :math:`v_i` is the velocity of the  :math:`i^\mathrm{th}` particle of the rigid body.
- :math:`\omega` is the angular velocity of the rigid body.
- :math:`I` is the inertia tensor, with components :math:`I_{i,j}`.
- :math:`L` is the angular momentum, a vector.
- :math:`T` is the torque on the rigid body. Elsewhere in this documentation, :math:`T` is
  the temperature. The context is sufficient to avoid ambiguity.
- :math:`t` is the time.
- :math:`h` is the integrator time step.
- :math:`\dot \bullet` is the time derivative of :math:`\bullet`.
- :math:`\wedge` is the cross product.
- :math:`\cdot` is the dot product.
- :math:`\hat x` is the unit vector parallel to :math:`x`.
- :math:`q`, :math:`q_1`, etc. are quaternions.
- :math:`q^\ast` is the conjugate of  :math:`q`.
- :math:`|q|` is the norm of :math:`q`.
- The superscript :math:`^B` is for quantities defined in the frame of reference of the
  rigid body. Other quantities are in the laboratory reference frame.


Quaternions are defined as the sum of a scalar part :math:`s` and a 3-dimensional vector part :math:`w`
as

.. math::
   q = s + w

The multiplication of the elementary vectors :math:`w_i=(1, 0, 0)`, :math:`w_j=(0, 1, 0)` and
:math:`w_k=(0,0,1)` is defined as

- :math:`w_i w_j = w_k`
- :math:`w_j w_k = w_i`
- :math:`w_k w_i = w_j`
- :math:`w_j w_i = -w_k`
- :math:`w_k w_j = -w_i`
- :math:`w_i w_k = -w_j`
- :math:`w_i^2 = w_j^2 = w_k^2 = -1`


From these definitions, one can derive the product of quaternions:

.. math::
   q_1 q_2 = \left(s_1 s_2 - w_1\cdot w_2\right) + \left(q_1 w_2 + q_2 w_1 + w_1 \wedge w_2\right)

where the first parenthesized term is the scalar part and the second one the vector
part. The addition of quaternions is done separately for the scalar and vector parts,
following conventional algebra. The conjugate of a quaternion is :math:`q^\ast = s - w`, its norm
is :math:`|q|=\sqrt{s^2 + w\cdot w}` and its inverse is :math:`q^{-1} = q^\ast / |q|^2`.

The utility of quaternions arises from their ability to encode any solid rotation. The
rotation operator about an axis :math:`n` with an angle :math:`\theta` is given by :math:`x' = q x q^\ast`
where :math:`q` is a unit quaternion with scalar component :math:`\cos\theta/2` and vector component
:math:`\sin\theta/2 \ \hat n`.

All quaternion operations in RMPCDMD are performed with the Fortran module
`fortran_quaternion <https://github.com/pdebuyl/fortran_quaternion>`_.

Definition of rigid-body quantities
-----------------------------------

The position in the laboratory frame of a member of the rigid body is

.. math::
   r_i(t) = q(t) r_i^B q^\ast(t)

The equation of motion for the angular momentum in the laboratory frame is

.. math::
    \dot L = T ,

and in the body frame

.. math::
    \dot L^B = T^B - \omega^B \wedge L^B .

The equation of motion for the quaternion is

.. math::
   \dot q = \frac{1}{2} \omega q


Algorithm
---------

Here, we follow section IV.A of Ref. :cite:`rozmanov_rotational_2010`.


+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| Angular momentum                                           | Torque	                             | Quaternion                                                       | Angular velocity                                          |
+============================================================+=======================================+==================================================================+===========================================================+
| :math:`L(h/2) = L(0) + h T(0)/2`                           |                                       |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| :math:`L^B(0) = q^\ast(0) L(0) q(0)`                       | :math:`T^B(0) = q^\ast(0) T(0) q(0)`  |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       |                                                                  | :math:`\omega^B(0) = L^B(0) / I^B`                        |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| :math:`\dot L^B(0) = T^B(0) - \omega^B(0)\wedge L^B(0)`    |                                       |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| :math:`L^B(h/2) = L^B(0) + h \dot L^B(0)/2`                |                                       |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       |                                                                  | :math:`\omega^B(h/2) = L^B(h/2)/I^B`                      |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       | :math:`^0\dot q(h/2) = q(0) \omega^B(h/2)/2`                     |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       | :math:`^0q(h/2) = q(0) + h ~^0\dot q(h/2)/2`                     |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| Start of the iterative procedure                                                                                                                                                                                                  |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| :math:`^{k+1}L^B(h/2) = ~^kq^\ast(h/2) L(h/2) ~^kq(h/2)`   |                                       |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       |                                                                  | :math:`^{k+1}\omega^B(h/2) = ~^{k+1}L^B(h/2) / I^B`       |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       | :math:`^{k+1}\dot q(h/2) = ~^k q(h/2) ~^{k+1}\omega^B(h/2)/2`    |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       | :math:`^{k+1}q(h/2) = q(0) + h~ ^{k+1}\dot q(h/2)/2`             |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|  End of the iterative procedure                                                                                                                                                                                                   |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            |                                       | :math:`q(h) = q(0) + h \dot q(h/2)`                              |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| Update positions using :math:`q(h)`                                                                                                                                                                                               |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| Compute the updated forces and torques in the lab frame                                                                                                                                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
|                                                            | :math:`T^B(h) = q^\ast(h)T(h)q(h)`    |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| :math:`L(h) = L(h/2) + h T(h)/2`                           |                                       |                                                                  |                                                           |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+
| Update the velocities using :math:`\omega^B(h) = L^B(h)/I^B`.                                                                                                                                                                     |
+------------------------------------------------------------+---------------------------------------+------------------------------------------------------------------+-----------------------------------------------------------+

The steps until the update of :math:`q(h)` and of the positions is the first step of the
velocity Verlet algorithm. These steps are implemented in :doxytag:`rigid_body_vv1`.

The update of :math:`L(h)` and of the velocities form the second part of the velocity Verlet
algorithm and are implemented in :doxytag:`rigid_body_vv2`.
