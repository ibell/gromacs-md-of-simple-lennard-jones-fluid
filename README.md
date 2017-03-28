# Building gromacs from sources


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
