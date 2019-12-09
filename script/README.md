# Scripts

This directory contains:

- The Material_Model folder. Within it are two matlab scripts to fit either a linear or a neo Hookean material model. Further explanation on how to use the scripts are contained in their header. The two mat files contain the test data from the mechanical characterization tests.

- batch_cellogram.m is a template matlab script to run cellogram in batch mode, which is substantially faster. Further instructions are contained within and can be taken from the cellogram documentation.

- read_json.m is a matlab script that takes the json formattted output from cellogram and converts it into a matlab structure for convenience of downstream processing.
