.. _matlab:

=============
 MATLAB code
=============

This MATLAB code is documented with Sphinx
using the `matlabdomain extension <https://github.com/sphinx-contrib/matlabdomain/blob/master/README.rst>`_.


src
========

.. mat:automodule:: src

:mod:`src` module contains the following source code files:
    
.. mat:autofunction:: src.main


src/integrators
===============

.. mat:autoscript:: src.integrators.LieEuler

.. mat:autoscript:: src.integrators.RKMK2Heun

.. mat:autoscript:: src.integrators.RKMK3

.. mat:autoscript:: src.integrators.RKMK4

.. mat:autoscript:: src.integrators.CFree4


src/Lie_group_functions
===============

.. mat:autoscript:: src.Lie_group_functions.hat

.. mat:autoscript:: src.Lie_group_functions.invhat

.. mat:autoscript:: src.Lie_group_functions.expSO3

.. mat:autoscript:: src.Lie_group_functions.expSE3

.. mat:autoscript:: src.Lie_group_functions.dexpinvSO3

.. mat:autoscript:: src.Lie_group_functions.dexpinvSE3


src/pre_post_processing
=======================

.. mat:autoscript:: src.pre_post_processing.preprocess

.. mat:autoscript:: src.pre_post_processing.postprocess

.. mat:autoscript:: src.pre_post_processing.plots


src/control_functions
=======================

.. mat:autoscript:: src.control_functions.control

.. mat:autoscript:: src.control_functions.getUpar

.. mat:autoscript:: src.control_functions.getUperp

.. mat:autoscript:: src.control_functions.deriv1

.. mat:autoscript:: src.control_functions.deriv2

.. mat:autoscript:: src.control_functions.deriv3

.. mat:autoscript:: src.control_functions.deriv4

.. mat:autoscript:: src.control_functions.derivhatq2w

.. mat:autoscript:: src.control_functions.deriv2hatq2w

.. mat:autoscript:: src.control_functions.deriv3hatq2w

.. mat:autoscript:: src.control_functions.deriv4hatq2w