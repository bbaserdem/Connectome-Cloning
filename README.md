# Connectome Cloning

This document countains  details regarding the connectome cloning project.

___
# Tasks

Here are some finalizing tasks.

* Flickering synapses (15% of all available synapses are present)

___
# Discussion

For flickering synapses, since both swapping and flipping is direction agnostic,
the flicker map must be symmetric.

___
# Workflow

The function **connclone** is the latest version of the algorithm.
That function can be checked to see the latest format to enter data.
All scripts in this directory are wrappers to record and save data.
Graph generating scripts are put in *visuals/generators*.

___
# Variables

Free letters remaining: GP

These are variables used to keep track of things.
Originally on **connclone**, these might differ between versions.
Will add version spoecific ones here as well.

## Inputs
* N **100**: Input connectome size
* D **0.1**: Connectome density
* E **10** : Cost parameter
* H **0,0**: Temperature scaling.
* *U* **0.3**: Active synapse percentage *(flickering only)*.

## Outputs
* O: 0 if reconstruction was unsuccessful, convergence time if successful.
* B: number of barcode pairs in the connectome

## Runtime
* W: Used for batch generating random numbers.
* C: The original connectome
* V: Reconstructed connectome
* T: Max simulation runtime
* R: Move probabilities

*Flicker*: for changing activation of synapses

* P: Number of times the synapse activation should change.

## State
* K: List of barcodes on bp ends (pre,pos) 
* S: List that contains bp locations. (pre,pos)

## Tracking
For hamiltonian;

* M: Table keeping track of barcodes per cell (barc,cell)
* F: Total barcodes per cell (cell)

For cell selection

* A: List of bp ids in each cell (cell,...)
* J: Total # of barcode pairs per cell (cell)
* Q: Location of bp in A (cell,bp-id)

For swap candidate selection

* Y: List of cells with >1 bp in them (...)
* Z: Total # of cells with >1 bp in them
* X: Location of cells in A (cell)

*Flicker*: for jump move selection

* I: List of synapses to cells.
* U: Total number of synapses per cell