# ParaGnosis

## Installation

A requirement is python3 and setuptool. These can be installed with (on debian based linux):

    > sudo apt-get install python3
    > sudo apt-get install python3-setuptools

The gears cli package can be installed with:

    > sudo python3 setup.py install

Or if you prefer to specify your own installation path <prefix>:

    > export PYTHONPATH=<prefix>:$PYTHONPATH
    > python3 setup.py install --prefix=<prefix>

## Running

Once installed you can launch a terminal and type the following to get a help message:

    > gears

You can also run the gears cli locally, without installation:

    > <path>/<to>/<gears-cli>/bin/gears


# compile asia.net in ../data/net/asia.net with tdmg and ACE

    > ./test/run.py --test compilation --network asia --bdd tdmg ace

# exhaustive marginal inference comparison with asia.net, with tdmg and ACE

    > ./test/run.py --test inference --network asia --bdd tdmg ace

# exhaustive marginal inference comparison with verification, with asia.net, with tdmg and ACE

    > ./test/run.py --test inference --verify --network asia --bdd tdmg ace


