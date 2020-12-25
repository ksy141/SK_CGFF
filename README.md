# SK_CGFF
Siyoung Kim's CG model

Run in gromacs via
gmx grompp -f step7.mdp -c yourstructure.gro -o step7.tpr
gmx mdrun -deffnm step7 -table toppar/table.xvg


Make topol.top by yourself (Write how many PL and TG molecules are in your system.)


toppar/DOPC.itp and toppar/TRIO.itp describe PL and TG molecules, respectively, defining the atom types, bonds, and angles. (No need to modify.)


toppar/cg.itp defines the non-bonded potential for each pair.
cd toppar
python make_potential.py
will create cg.itp.


The initial structures for a PL bilayer or PL+TG bilayer can be made with the attached scripts.
python build_PL.py or python build_PLTG.py.
Requires the MDAnalysis package.
