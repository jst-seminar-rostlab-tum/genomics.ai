Genecruncher Documentation
============================

Genecruncher helps you visualize all of your single-cell sequencing data in a fast and easy way using neural networks.
The neural networking is done using the scArches package. More information about scArches can be found `here <https://scarches.readthedocs.io/en/latest/>`_.

First steps
---------

First login or sign up using your academic affiliation.

.. image:: _static/homepage.png

Under the Gene Mapper all of your mappings are displayed. Click on the plus button to create a new one.

.. image:: _static/mappings_dashboard.png

On this page you should pick an Atlas and a model.
Please note that some models are not compatible with some atlases. If so, they will be disabled. 
Here the Atlas and model type are set, scVI for the Pancreas for example.

.. image:: _static/new_mapping_1.png

On the last page the query datasets are selected. 
Here you are able to review them and go back via the back button or by clicking on the step indicator should you change your mind.
Right below that the consequent requirements are listed as a result of your choices.
These requirements should be kept in mind if you wish to upload your own dataset, which could be done by drag-and-drop or by clicking on the upload field.
Alternatively you could select one of the existing ones.

.. image:: _static/new_mapping_2.png

Click on "Create Mapping", give it an appropriate name and submit.
You will then be forwarded back to the homepage where you could eventually review the result once it is done processing.
Please note that it takes some time to process.

A color-coded build status will appear with a progress bar. As soon as it's ready you can click on "See Results".

Here the mapping in the center can be moved and resized via mouse move and scroll, the color-coding modes and selections are toggleable via the Categories tab on the left and on the right the displayed data is visualized in a barchart to display differences between the different batch types and other accumulative values. On the top there are a few buttons to toggle the visability of reference and query data and more.

.. image:: _static/visualizer.png


Links To The Detailed Information Pages
---------

.. toctree::
   :maxdepth: 1

   index.rst
   visualization/index.rst