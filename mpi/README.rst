---
mpi
---

mpi contains several modules build to interface with HPC systems using the installed MPI library (through mpi4py). This comes in two pieces: parallelization of frame processing and accumulation.

Frame IO Parallelization
========================

The MPIFileIO sub-module contains a few modules that will add frame-level parallelism to a G3Pipeline by reading files within an observation in parallel and (optionally, with MPIIODistributor or MPIFrameParallelizer) distributes the data frame-by-frame across any number of processes with an M:N IO:CPU node distribution.

This processing will make sure all metadata frames are seen by all processes in the same order with respect to both each other and all data frames and are thus strongly ordered across all processes. Data frames (scans, typically) are ordered on a node, but are not ordered between nodes and thus should be considered weakly ordered. For a frame sequence MN12345PQ67, where letters denote metadata frames (calibration, etc.) and numbers data frames, a possible ordering seen by the pipelines on two processing nodes would be MN134PQ7 and MN25PQ6, respectively. A consequence of these ordering rules is that you *must not use any modules that depend on frame ordering and continuity*. Examples of such would be modules that buffer multiple scans together to get longer-time-period data; by their nature, these only work for observation-level parallelism, rather than frame-level or file-level.

Parallelization in this sense can be achieved by replacement of G3Reader in a normal G3Pipeline with MPIIODistributor.

Frame Accumulation
==================

At the end of a sequence of parallel pipelines, you in general want to join the data from all the processes in the communicator back together again. Typical cases are stitching processed scans back together again, for example to feed to a maximum-likelihood map-maker or some other algorithm that needs very large quantities of distributed data, or doing a parallel reduction, for example by coadding map frames produced by N parallel map-makers.

The first of these cases (restitching an observation processed on many parallel processes) is more complicated than the second, so this library provides a helper module (MPIAccumulator) that does it for you, placing the result in a member variable inspected after the pipelines end.

The second case is easier, but potentially common enough that a module should be added in future.

Interface to TOAST
==================

There is an experimental module (TOASTFiller) that uses MPIAccumulator to fill frame contents after processing in many parallel pipelines into a TOAST TOD object in order to use TOAST's set of parallel timestream algorithms and map-making. Basic functionality is present, but it should be considered in-development at present.
