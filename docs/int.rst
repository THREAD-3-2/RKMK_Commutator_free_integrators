.. _int:

==========================================
RKMK and Commutator free Lie group methods
==========================================

The underlying idea of RKMK methods is to express a vector field :math:`F\in\mathfrak{X}(\mathcal{M})` as :math:`F\vert_m = \Psi_*(f(m))\vert_m` , where :math:`\Psi_*` is the infinitesimal generator of :math:`\Psi`, a transitive action on :math:`\mathcal{M}`, and :math:`f:\mathcal{M}\rightarrow\mathfrak{g}`. This allows us to transform the problem from the manifold :math:`\mathcal{M}` to the Lie algebra :math:`\mathfrak{g}`, on which we can perform a time step integration. We then map the result back to :math:`\mathcal{M}`, and repeat this up to the final integration time.  More explicitly, let :math:`h_n` be the size of the :math:`n-th` time step, we then update :math:`y_n\in\mathcal{M}` to :math:`y_{n+1}` by

.. math::
    :name: eq:1
    
    \begin{align}
        \begin{cases}
        \sigma(0) = 0\in\mathfrak{g},\\
        \dot{\sigma}(t) = \textrm{dexp}_{\sigma(t)}^{-1}\circ f\circ \Psi (\exp(\sigma(t)),y_n)\in T_{\sigma(t)}\mathfrak{g}, \\
        y_{n+1} = \Psi(\exp(\sigma_1),y_n)\in \mathcal{M},
        \end{cases}
    \end{align}

where :math:`\sigma_1\approx \sigma(h_n)\in\mathfrak{g}` is computed with a Runge-Kutta method. To be finished...
