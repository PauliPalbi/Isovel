# Isovel

Isovel is a python package that draws isovelocity curves over gas images of [Disk Substructures at High Angular Resolution Project (DSHARP) ](https://almascience.nrao.edu/almadata/lp/DSHARP/)conducted with the Atacama Large Milimeter/Submilimeter Array (ALMA).

## Getting Started
The most simple way is to use  `pip`,

```
pip install Isovel
```
Or by cloning this repository

```bash
# Clone this repository
git clone https://github.com/PauliPalbi/Isovel

# Go into the repository
cd Isovel

# Install required modules
pip install -r requirements.txt
```

## Running a Test
Start by using an fiducial image (_CO.fits) from the [DSHARP Data Release webpage](https://almascience.nrao.edu/almadata/lp/DSHARP/images/).

In the case where the shape of the gas of the disk is unknown, install [`eddy`](https://github.com/richteague/eddy)

```
pip install astro-eddy
```

and [bettermoments](https://github.com/richteague/bettermoments) 

```
pip install bettermoments

bettermoments path/to/cube.fits
```


## Acknowledgments

* This project is part of the [2020 Code/Astro workshop](https://semaphorep.github.io/codeastro/) to develope and release an open-source astromony package.