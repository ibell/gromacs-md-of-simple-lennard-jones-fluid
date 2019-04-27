# Building gromacs from sources

All informations here: http://manual.gromacs.org/documentation/2016.3/install-guide/index.html
```sh
git clone https://github.com/gromacs/gromacs
cd gromacs
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON
make -j
make check
sudo make -j install
source /usr/local/gromacs/bin/GMXRC
```
You should add the "source" to your bashrc or zshrc or ... that gives you access to gmx through


# Generate a lattice of 9000 methane (single site) molecules

```sh
./generate_lattice.py --density=0.668 --atom-name=Me --lattice=cubic --molecular-weight=16.0425  9000 
```

Note that more than 10k atoms seems to crash the pdb (surely pb with columns requirements of pdf).
The script has been found at https://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/11/SimulatingliquidArgonwithGromacs.pdf

```sh
gmx grompp -f em.mdp -c methane_start.pdb -p methane.top -o methane_em.tpr # generates the tpr file to be used by mdrun
gmx mdrun -s methane_em.tpr -v -c methane_em.gro # generates the trajectory file
gmx grompp -f npt_md.mdp -c methane_em.gro -p methane.top -o methane_npt.tpr
gmx mdrun -s methane_npt.tpr -v -c methane_npt.gro
```

# Useful ressources

* https://encyclopedia.airliquide.com/fr/methane
* http://manual.gromacs.org/programs/gmx-grompp.html
* http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/
* https://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/11/SimulatingliquidArgonwithGromacs.pdf
* http://www3.mpibpc.mpg.de/groups/de_groot/compbio/p1/index.html
* http://manual.gromacs.org/documentation/2016.3/install-guide/index.html
* https://github.com/gromacs/gromacs
