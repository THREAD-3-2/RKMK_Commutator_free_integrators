.. _int:

==========================================
RKMK and Commutator free Lie group methods
==========================================

Lie group integrators solve differential equations whose solution evolve on a manifold :math:`\mathcal{M}`, i.e. the solution is a curve :math:`y(t)\in\mathcal{M}` whose tangent at any point coincides with a vector field :math:`F\in\mathcal{X}(\mathcal{M})` and passes through a designated initial value :math:`y_0` at :math:`t=t_0`:

.. math::

    \begin{align}
        \dot{y}(t) = F|_{y(t)},\qquad y(t_0)=y_0.    
    \end{align}

The underlying idea of Runge-Kutta-Munthe-Kaas (RKMK) methods is to express a vector field :math:`F\in\mathfrak{X}(\mathcal{M})` as :math:`F\vert_m = \Psi_*(f(m))\vert_m` , where :math:`\Psi_*` is the infinitesimal generator of :math:`\Psi`, a transitive action on :math:`\mathcal{M}`, and :math:`f:\mathcal{M}\rightarrow\mathfrak{g}`. This allows us to transform the problem from the manifold :math:`\mathcal{M}` to the Lie algebra :math:`\mathfrak{g}`, on which we can perform a time step integration with a Runge-Kutta method. We then map the result back to :math:`\mathcal{M}`, and repeat this up to the final integration time.  More explicitly, let :math:`h_n` be the size of the :math:`n-th` time step, we then update :math:`y_n\in\mathcal{M}` to :math:`y_{n+1}` by

.. math::
    :name: eq:1
    
    \begin{align}
        \begin{cases}
        \sigma(0) = 0\in\mathfrak{g},\\
        \dot{\sigma}(t) = \textrm{dexp}_{\sigma(t)}^{-1}\circ f\circ \Psi (\exp(\sigma(t)),y_n)\in T_{\sigma(t)}\mathfrak{g}, \\
        y_{n+1} = \Psi(\exp(\sigma_1),y_n)\in \mathcal{M},
        \end{cases}
    \end{align}

where :math:`\sigma_1\approx \sigma(h_n)\in\mathfrak{g}` is computed with a Runge-Kutta method. 


The transformed differential equation for :math:`\sigma(t)` makes use of the derivative of the exponential mapping. The map :math:`v\mapsto\textrm{dexp}_u(v)` is linear and invertible when :math:`u` belongs to some sufficiently small neighborhood of :math:`0\in\mathfrak{g}`. It has an expansion in nested Lie brackets and, using the operator :math:`\textrm{ad}_u(v)=[u,v]` and its powers :math:`\textrm{ad}_u^2 v=[u,[u,v]]` etc, one can write

.. math::
    :name: eq:2
    
    \begin{align}
        \textrm{dexp}_u(v) = \left.\frac{e^z-1}{z}\right|_{z=\textrm{ad}_u}(v) = v + \frac12[u,v] + \frac16[u,[u,v]] + \cdots.
    \end{align}

The inverse is

.. math::
    :name: eq:3
    
    \begin{align}
        \textrm{dexp}_u^{-1}(v) =\left.\frac{z}{e^z-1}\right|_{z=\textrm{ad}_u}(v)= v -\frac12[u,v] + \frac1{12}[u,[u,v]]+\cdots.
    \end{align}

To evaluate :math:`\textrm{dexp}_u^{-1}(v)` one can either truncate the series :ref:`(3) <eq:2>`, or compute its exact expression for the particular Lie algebra under consideration.

The second class of Lie group integrators to be considered here are the commutator-free methods, named this way to emphasize the contrast to RKMK schemes which usually include commutators in the method format. These schemes include the Crouch-Grossman methods and have the format

.. math::
    
    \begin{align}
        Y_{n,r} &= \exp\left(h\sum_{k}\alpha_{r,J}^k f_{n,k}\right)\cdots \exp\left(h\sum_{k}\alpha_{r,1}^k f_{n,k}\right)\cdot y_n \\
        f_{n,r} &= f(Y_{n,r}) \\[1mm]
        y_{n+1} &= \exp\left(h\sum_k \beta_J^k f_{n,k}\right)\cdots \exp\left(h\sum_k \beta_1^k f_{n,k}\right)\cdot y_n
    \end{align}

Here the Runge-Kutta coefficients :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are related to a classical Runge-Kutta scheme with coefficients :math:`a_r^k`, :math:`b_r` in that :math:`a_r^k=\sum_j \alpha_{r,j}^k` and :math:`b_r=\sum_j \beta_{j}^r`. The :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are usually chosen to obtain computationally inexpensive schemes with the highest possible order of convergence. The computational complexity of the above schemes depends on the cost of computing an exponential as well as of evaluating the vector field. Therefore it makes sense to keep the number of exponentials :math:`J` in each stage as low as possible, and possibly also the number of stages :math:`s`.

We refer to `(Celledoni, Çokaj, Leone, Murari and Owren, 2021) <https://doi.org/10.1080/00207160.2021.1966772>`_ and references therein for further details.




