# BXSF and FRMSF File Handling

This repository contains Python functions and classes for handling data from BXSF and FRMSF files commonly used in materials science and condensed matter physics. The provided scripts enable you to read, manipulate, and write data from these file formats.

## Table of Contents

- [Introduction](#introduction)
- [File Handling](#bxsf-file-handling)

## Introduction

BXSF and FRMSF files are commonly used for storing data related to electronic band structures, Fermi surfaces, and other materials properties. This repository provides Python functions and classes to work with these file formats.

## BXSF and FRMSF File Handling

The `fermio` module in this repository contains `bxsf` and `frmsf` for reading and writing BXSF/frmsf files. You can use this class to read data from files, manipulate the data, and write the modified data back to either file. Here are some of the features:

- Read BXSF files and access the lattice vectors, energy bands, Fermi energy, and k-point data.
- Write data to a new BXSF/FRMSF file.

```python
from fermio import bxsf

# Read a BXSF file
bxsf_obj = bxsf.from_file('sample.bxsf')

# Access data and manipulate it
k_vals = bxsf_obj.k_list

# Write to a new BXSF file
bxsf_obj.write('new_sample.bxsf')
