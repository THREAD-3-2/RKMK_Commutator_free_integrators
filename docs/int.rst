.. _int:

===========
Algorithms
===========

.. _rkmk_int:

Runge-Kutta-Munthe-Kaas methods
-------------------------------

Lie group integrators solve differential equations whose solution evolve on a 
manifold :math:`\mathcal{M}`, i.e. the solution is a curve :math:`y(t)\in\mathcal{M}` 
whose tangent at any point coincides with a vector field :math:`F\in\mathcal{X}(\mathcal{M})` 
and passes through a designated initial value :math:`y_0` at :math:`t=t_0`:

.. math::
    :label: vecfield

    \begin{align}
        \dot{y}(t) = F\left(y(t)\right),\qquad y(t_0)=y_0.    
    \end{align}

Notice that we restrict the discussion to the case of autonomous vector field (explicit time dependence 
could easily be included). Let :math:`G` be a Lie group acting transitively on :math:`\mathcal{M}` via the 
group action :math:`\psi:G \times \mathcal{M} \rightarrow \mathcal{M}`, so 
that :math:`\mathcal{M}` is a homogeneous manifold. The underlying idea of Runge-Kutta-Munthe-Kaas 
(RKMK) methods is to express a vector field :math:`F\in\mathfrak{X}(\mathcal{M})` as 
:math:`F\vert_m = \psi_*(f(m))\vert_m` , where :math:`\psi_*` is the infinitesimal generator 
of :math:`\psi` and :math:`f:\mathcal{M}\rightarrow\mathfrak{g}`. This allows us to transform 
the problem from the manifold :math:`\mathcal{M}` to the Lie algebra :math:`\mathfrak{g}` 
of :math:`G`, on which we can perform a time step integration with a Runge-Kutta method. We 
then map the result back to :math:`\mathcal{M}`, and repeat this up to the final integration time. 
More explicitly, let :math:`h_n` be the size of the :math:`n-th` time step, we then update 
:math:`y_n\in\mathcal{M}` to :math:`y_{n+1}` by

.. math::
    :label: int2
    
    \begin{align}
        \begin{cases}
        \sigma(0) = 0\in\mathfrak{g},\\
        \dot{\sigma}(t) = \textrm{dexp}_{\sigma(t)}^{-1}\circ f\circ \psi \left(\exp(\sigma(t)),y_n\right)\in T_{\sigma(t)}\mathfrak{g}, \\
        y_{n+1} = \psi(\exp(\sigma_1),y_n)\in \mathcal{M},
        \end{cases}
    \end{align}

where  :math:`\textrm{exp}:\mathfrak{g}\rightarrow G` is the exponential map, 
and :math:`\sigma_1\approx \sigma(h_n)\in\mathfrak{g}` is computed with a Runge-Kutta method. 


The transformed differential equation for :math:`\sigma(t)` makes use of the derivative of 
the exponential mapping. The map :math:`v\mapsto\textrm{dexp}_u(v)` is linear and invertible 
when :math:`u` belongs to some sufficiently small neighborhood of :math:`0\in\mathfrak{g}`. It 
has an expansion in nested Lie brackets and, using the operator :math:`\textrm{ad}_u(v)=[u,v]` 
and its powers :math:`\textrm{ad}_u^2 v=[u,[u,v]]` etc, one can write

.. math::
    :label: int3
    
    \begin{align}
        \textrm{dexp}_u(v) = \left.\frac{e^z-1}{z}\right|_{z=\textrm{ad}_u}(v) = v + \frac12[u,v] + \frac16[u,[u,v]] + \cdots.
    \end{align}

The inverse is

.. math::
    :label: int4
    
    \begin{align}
        \textrm{dexp}_u^{-1}(v) =\left.\frac{z}{e^z-1}\right|_{z=\textrm{ad}_u}(v)= v -\frac12[u,v] + \frac1{12}[u,[u,v]]+\cdots.
    \end{align}

To evaluate :math:`\textrm{dexp}_u^{-1}(v)` one can either truncate the series :ref:`(3) <eq:2>`, 
or compute its exact expression for the particular Lie algebra under consideration. The exponential 
maps on the Lie groups SO(3) and SE(3) are implemented in 
`expSO3 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/lie_group_functions/expSO3.m>`_ 
and `expSE3 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/lie_group_functions/expSE3.m>`_. 
The exact expressions for the inverse of the derivative of 
the exponential map on SO(3) and SE(3) are implemented 
in `dexpinvSO3 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/lie_group_functions/dexpinvSO3.m>`_ 
and `dexpinvSE3 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/lie_group_functions/dexpinvSE3.m>`_.

.. _rkmk_ex_int:

Examples
^^^^^^^^
Let us consider an s-stage Runge-Kutta (RK) method for the initial value problem :eq:`vecfield`:

.. math::
    :label: int5

    \begin{align}
    y_{n+1}=y_n+h \sum_{i=1}^s b_i k_i, \quad k_i=F\left(y_n+h \sum_{j=1}^s a_{i j} k_j \right), \quad i=1, \ldots, s,
    \end{align}

where :math:`b_i,\,a_{ij}\, (i,\,j=1,\dots\,s)` are real numbers called, respectively, the weights and coefficients of 
the method, and :math:`c_i=\sum_{j=1}^s a_{ij}` are called the nodes or abscissae. These constants define a specific RK method and can 
be collected in the following table, known as Butcher's tableau:

.. math::
    :label: int6

    \begin{align}
    \begin{array}{c|cccc}
    c_1 & a_{11} & a_{12} & \ldots & a_{1 s} \\
    c_2 & a_{21} & a_{22} & \ldots & a_{2 s} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    c_s & a_{s 1} & a_{s 2} & \ldots & a_{s s} \\
    \hline & b_1 & b_2 & \ldots & b_s
    \end{array}
    \end{align}

From equation :eq:`int2` it follows that one step of the resulting Runge–Kutta–Munthe-Kaas method writes 

.. math::
    :label: int7

    \begin{align}
    &y_1=\exp \left(h \sum_{i=1}^s b_i k_i\right) \cdot y_0,\\
    &k_i=\operatorname{dexp}^{-1}_{h \sum_j a_{i j} k_j} f\left(\exp \left(h \sum_j a_{i j} k_j\right) \cdot y_0\right), \quad i=1, \ldots, s,
    \end{align}

where we denote the group action by ":math:`\cdot`" for ease of notation. 

The simplest Lie group integrator is the Lie-Euler method, based on the classical explicit Euler method, a first-order method with Butcher's tableau given by

.. math::

    \begin{align}
    \begin{array}{c|c}
    0 & 0 \\
    \hline & 1
    \end{array}
    \end{align}

The resulting Lie-Euler method can be written as :math:`y_{n+1}=\exp \left(h F(y_n)\right) y_n` and is implemenmted in 
`LieEuler <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/integrators/LieEuler.m>`_.

An improvement to the Lie-Euler method is the second-order RKMK method based on the tableau of the Heun's method:

.. math::

    \begin{align}
    \begin{array}{c|cc}
    0 & 0 & 0 \\
    1 & 1 & 0 \\
    \hline & 1/2 & 1/2
    \end{array}
    \end{align}

The resulting RKMK integrator is implemented in `RKMK2Heun <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/integrators/RKMK2Heun.m>`_

The following Butcher's tables provide the coefficients for two classical methods of order three (on the left) and order four (on the right):

.. math::
    :label: 4ord

    \begin{align}
    \begin{array}{c|ccc}
    0 & 0 & 0 & 0 \\
    {1/2} & {1/2} & 0 & 0 \\
    1 & -1 & 2 & 0 \\
    \hline & {1/6} & {2/3} & {1/6}
    \end{array} 
    \qquad \qquad \quad
    \begin{array}{c|cccc}
    0 & 0 & 0 & 0 & 0 \\
    {1/2} & {1/2} & 0 & 0 & 0 \\
    {1/2} & 0 & {1/2} & 0 & 0 \\
    1 & 0 & 0 & 1 & 0 \\
    \hline & {1/6} & {1/3} & {1/3} & {1/6}
    \end{array}
    \end{align}

The corresponding RKMK integrators are implemented in `RKMK3 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/integrators/RKMK3.m>`_ and
`RKMK4 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/integrators/RKMK4.m>`_.




.. _cfree_int:

Commutator-free methods
-----------------------

The second class of Lie group integrators to be considered here are the commutator-free methods, 
named this way to emphasize the contrast to RKMK schemes which usually include commutators in 
the method format. These schemes include the Crouch-Grossman methods and have the format

.. math::
    
    \begin{align}
        Y_{n,r} &= \exp\Big(h\sum_{k}\alpha_{r,J}^k f_{n,k}\Big)\cdots \exp\Big(h\sum_{k}\alpha_{r,1}^k f_{n,k}\Big) \cdot y_n\\
        f_{n,r} &= f(Y_{n,r}) \\[1mm]
        y_{n+1} &= \exp\Big(h\sum_k \beta_J^k f_{n,k}\Big)\cdots \exp\Big(h\sum_k \beta_1^k f_{n,k}\Big) \cdot y_n
    \end{align}

where ":math:`\cdot`" denotes the group action. Here 
the Runge-Kutta coefficients :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are related 
to a classical Runge-Kutta scheme with coefficients :math:`a_r^k`, :math:`b_r` in 
that :math:`a_r^k=\sum_j \alpha_{r,j}^k` and :math:`b_r=\sum_j \beta_{j}^r`. 
The :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are usually chosen to obtain 
computationally inexpensive schemes with the highest possible order of convergence. 
The computational complexity of the above schemes depends on the cost of computing an 
exponential as well as of evaluating the vector field. Therefore it makes sense to 
keep the number of exponentials :math:`J` in each stage as low as possible, and 
possibly also the number of stages :math:`s`. 

The following example is a generalization of the classical fourth-order
Runge–Kutta method in :eq:`4ord` and is implemented in  `CFree4 <https://github.com/THREAD-3-2/RKMK_Commutator_free_integrators/blob/main/src/integrators/CFree4.m>`_:


.. math::

    \begin{aligned}
    &Y_1=y_0, \\
    &Y_2=\exp \left(\frac{1}{2} k_1\right) \cdot y_0, \\
    &Y_3=\exp \left(\frac{1}{2} k_2\right) \cdot y_0 \\
    &Y_4=\exp \left(k_3-\frac{1}{2} k_1\right) \cdot Y_2, \\
    &y_{\frac{1}{2}}=\exp \left(\frac{1}{12}\left(3 k_1+2 k_2+2 k_3-k_4\right)\right) \cdot y_0, \\
    &y_1=\exp \left(\frac{1}{12}\left(-k_1+2 k_2+2 k_3+3 k_4\right)\right) \cdot y_{\frac{1}{2}},
    \end{aligned}

with :math:`k_i=hf(Y_i)`. We notice that one exponential is saved in computing :math:`Y_4` by making use of :math:`Y_2`. This shows that sometimes it is possible 
to come up with tricks that allow to reuse exponentials from one stage to another, thereby lowering 
the computational cost of the scheme.

We refer to `(Celledoni, Çokaj, Leone, Murari and Owren, 2021) <https://doi.org/10.1080/00207160.2021.1966772>`_ 
and references cited therein for further details.




