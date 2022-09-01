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
    :name: eq:int1
    \begin{align}
        \dot{y}(t) = F|_{y(t)},\qquad y(t_0)=y_0.    
    \end{align}

Let :math:`G` be a Lie group acting transitively on :math:`\mathcal{M}` via the 
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
    :name: eq:int2
    
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
    :name: eq:int3
    
    \begin{align}
        \textrm{dexp}_u(v) = \left.\frac{e^z-1}{z}\right|_{z=\textrm{ad}_u}(v) = v + \frac12[u,v] + \frac16[u,[u,v]] + \cdots.
    \end{align}

The inverse is

.. math::
    :name: eq:int4
    
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

One step of 

.. math::
    :name: eq:int5

\begin{align}
&y_1=\exp \left(h \sum_{i=1}^s b_i k_i\right) \cdot y_0,\\
&k_i=\operatorname{dexp}_{h \sum_j^{-1} a_{i j} k_j} f\left(\exp \left(h \sum_j a_{i j} k_j\right) \cdot y_0\right), \quad i=1, \ldots, s .
\end{align}


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

where we denote the group action by ":math:`\cdot`" for ease of notation. Here 
the Runge-Kutta coefficients :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are related 
to a classical Runge-Kutta scheme with coefficients :math:`a_r^k`, :math:`b_r` in 
that :math:`a_r^k=\sum_j \alpha_{r,j}^k` and :math:`b_r=\sum_j \beta_{j}^r`. 
The :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are usually chosen to obtain 
computationally inexpensive schemes with the highest possible order of convergence. 
The computational complexity of the above schemes depends on the cost of computing an 
exponential as well as of evaluating the vector field. Therefore it makes sense to 
keep the number of exponentials :math:`J` in each stage as low as possible, and 
possibly also the number of stages :math:`s`.

We refer to `(Celledoni, Çokaj, Leone, Murari and Owren, 2021) <https://doi.org/10.1080/00207160.2021.1966772>`_ 
and references therein for further details.




