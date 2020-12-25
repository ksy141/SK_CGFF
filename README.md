# SK_CGFF
Siyoung Kim's CG model

1. Run in gromacs via
``` bash
gmx grompp -f step7.mdp -c yourstructure.gro -o step7.tpr
gmx mdrun -deffnm step7 -table toppar/table.xvg
```

2. Make *topol.top* by yourself (Write how many PL and TG molecules are in your system.)

3. *toppar/DOPC.itp* and *toppar/TRIO.itp* describe PL and TG molecules, respectively, defining the atom types, bonds, and angles. (No need to modify.)

4. *toppar/cg.itp* defines the non-bonded potential for each pair. The below will make a new *cg.itp*.
``` bash
cd toppar
python make_potential.py
```

5. The initial structures for a PL bilayer or PL+TG bilayer can be made with the attached scripts (require MDAnalysis).
``` bash
python build_PL.py
python build_PLTG.py
```
