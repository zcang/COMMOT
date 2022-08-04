.. commot documentation master file, created by
   sphinx-quickstart on Sat Feb 20 12:08:49 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: commot
.. automodule:: commot
   :noindex:

API
==================================

Preprocessing: pp
-------------------

.. module:: commot.pp
.. currentmodule:: commot

.. autosummary::
   :toctree: .

   pp.infer_spatial_information
   pp.ligand_receptor_database
   pp.filter_lr_database


Tools: tl
-----------

.. module:: commot.tl
.. currentmodule:: commot

Communication inference
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.spatial_communication
   tl.communication_direction
   tl.cluster_communication
   tl.cluster_communication_spatial_permutation
   tl.cluster_position

Downstream analysis
~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.communication_deg_detection
   tl.communication_deg_clustering
   tl.communication_impact
   tl.group_cluster_communication
   tl.group_cell_communication
   tl.group_communication_direction
   tl.communication_spatial_autocorrelation


Plotting: pl
------------

.. module:: commot.pl
.. currentmodule:: commot

.. autosummary::
   :toctree: .

   pl.plot_cell_communication
   pl.plot_cluster_communication_network
   pl.plot_cluster_communication_dotplot
   pl.plot_communication_dependent_genes
   pl.plot_communication_impact
