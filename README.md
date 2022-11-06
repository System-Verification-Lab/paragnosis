# ParaGnosis

ParaGnosis is a weighted model counting toolset. Its implementation is based on [[1,2,3,4]](#4). We have also added a significant number of Bayesian networks to play with (under *./data/net*)

The toolset consists of the following:

  * `bn-to-cnf`: a c++ tool to create Conjunctive Normal Form (CNF) encodings from a Bayesian network.
  * `bnc`: a c/c++ **B**ayesian **N**etwork **C**ompiler for multiple target representations.
  * `bnmc`: a  c++ **B**ayesian **N**etwork **M**odel **C**ounter.
  * `pg`: a **P**ara**G**nosis user friendly jjinterface to the tools above, written in Python.

The currently supported target languages are:

  * Weighted Positive Binary Decision Diagrams (WPBDD)
  * Weighted Positive Multi-Valued Decision Diagrams (WPMDD)
  * Tree-driven Weighted Positive Multi-valued Decision Diagrams (TD-WPMDD)

## Installation (Ubuntu 18.04+)

Install requirements with `apt`:

    > sudo apt-get install -y libboost-all-dev python3 python-setuptools make cmake gcc g++ libgmp-dev libgsl-dev libreadline-dev make cmake

Install latest pip (the python package installer):

    > sudo python3 -m pip install --upgrade pip

To build all tools in the toolset, type:

    > make

Binaries will be installed in the `<path/to/source>/bin` directory, and the `pg` script will be available system wide.

## Reconfigure `pg`

In order to let the `pg` script know where the toolset is located, we can run `pg` commands with `pg --source-dir=<path/to/source> ...`, or adjust the following in `pg`'s configuration file `~/.pgrc`:

    location = <path/to/source>

To test if the installation is successful, we can give the following a try. Open a terminal at any location and type:

    > pg encode asia
    or:
    > pg --source-dir=<path/to/source> encode asia

This should produce encoding statistics for the *asia* network.

## Examples

For the following examples, we assume `location` has been set in `~/pgrc`. All available commands can be found through `pg --help`.

### Encoding and available networks

#### Show a list of available Bayesian networks

    > pg --list

        3nt
        4sp
        6hj
        6nt
        aggregate
        alarm
        ...

Any of the shown names can be used as input for the `pg` script.

#### Show encoding statistics

    > pg encode asia

        ...

        Variables       : 8
        Probabilities   : 36
        Deterministic   : 8
        Unsatisfiable   : 4
        Literals        : 16
        Clauses         : 52
        Literal/clauses : 2.23
        Clause sizes    : 1-3


### Compiling a network

#### Compile asia to a TD-WPMDD

    > pg compile asia

#### Compile asia to a WPBDD

    > pg compile asia --method wpbdd

#### Compare compilation between WPBDD, a WPMDD, and a TD-WPMDD

    > pg compile asia --method wpbdd mg tdmg

#### Compile a local Hugin file to

    > pg compile <path/to/file>/asia.net


#### Directly visualize (with evince) the WPBDD

    > pg compile asia --method wpbdd --dot

#### Directly visualize the TD-WPMDD

    > pg compile asia --method tdmg --dot

### Perform Inference

#### Run every possible marginalization on a network using TD-WPMDD (press ctrl-c to stop)

    > pg inference asia

#### Compute all posteriors for evidence `bronc = yes`, and `smoke = yes`.

    > pg inference asia --evidence='bronc=yes,smoke=yes'

#### Compute all posteriors of `lung` and `xray` for evidence `bronc = yes`, and `smoke = yes`

    > pg inference asia --evidence='bronc=yes,smoke=yes' --posteriors=lung,xray

        ...

        lung=yes: 0.100000
        lung=no: 0.900000
        xray=yes: 0.151705
        xray=no: 0.848295


## References

<a id="1">[1]</a>
G.H. Dal, A.W. Laarman, A. Hommerso and P.J.F. Lucas, ”*A Compositional Approach to Probabilistic Knowledge Compilation*”, in International Journal of Approximate Reasoning, vol 138:38-66, 2021.

<a id="2">[2]</a>
G.H. Dal, A.W. Laarman and P.J.F. Lucas, ”*Parallel Probabilistic Inference by Weighted Model Counting*”, in Proceeding of the International Conference on Probabilistic Graphical Models, PMLR, vol 72:97-108, 2018.

<a id="3">[3]</a>
G.H. Dal, S. Michels and P.J.F. Lucas, ”*Reducing the Cost of Probabilistic Knowledge Compilation*”, in Proceedings of Machine Learning Research, volume 73, pages 41-152, 2017.

<a id="4">[4]</a>
G.H. Dal and P.J.F. Lucas, ”*Weighted Positive Binary Decision Diagrams for Exact Probabilistic Inference*”, in Journal of Approximate Reasoning, volume 90, pages 411-432, 2017.
