# SalainDNA: A Pre-alignment Filter for DNA Read Mapping on a Multi-core Environment

## **Prerequisites**

To run the system, make sure that the following packages are installed:

- OpenMP
- GCC
- CMake

To install the prerequisites in Ubuntu, update first the Ubuntu repository indexes by running the command:

```bash
sudo apt update
```

Once the Ubuntu repository indexes are updated, run the following commands:

```bash
sudo apt install libomp-dev # for OpenMP
sudo apt install build-essentials # for GCC
sudo apt install cmake # for CMake
```



## Installation

To install, build the latest version found at the GitHub repository. In order to do this, run the following commands:

```bash
git clone https://github.com/sairakaye/salainDNA
cd salainDNA
cmake .
make
./salainDNA
```

After executing the last command, the name and steps included in the system as well as the set of arguments should show on the console.
