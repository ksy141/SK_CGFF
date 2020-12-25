import MDAnalysis as mda
import numpy as np

pbcz  = 150
N     = 100 * 2
Nx    = int(np.sqrt(N/2))
dx    = 8.5

sigma  = 5 #A
PL_natoms = 11

n_atoms = N*PL_natoms
n_residues = N
atom_resindex = np.repeat(np.arange(0, N), PL_natoms)

u = mda.Universe.empty(n_atoms = n_atoms,
                       n_residues = n_residues,
                       atom_resindex = atom_resindex,
                       trajectory = True)

u.add_TopologyAttr('name', ['NC3','PO4','PGL',
                            'C21', 'C22', 'C23', 'C24',
                            'C31', 'C32', 'C33', 'C34'] * N)

u.add_TopologyAttr('type', ['NC3','PO4','PGL',
                            'CT1','CT1','CT2','CT2',
                            'CT1','CT1','CT2','CT2'] * N)

u.add_TopologyAttr('resname', ['DOPC'] * N)
u.add_TopologyAttr('resid', np.arange(1, N+1))

bonds = []
for o in range(0, N*11, 11):
    bonds.append([o+0, o+1])
    bonds.append([o+1, o+2])
    bonds.append([o+2, o+3])
    bonds.append([o+2, o+7])
    bonds.append([o+3, o+4])
    bonds.append([o+4, o+5])
    bonds.append([o+5, o+6])
    bonds.append([o+7, o+8])
    bonds.append([o+8, o+9])
    bonds.append([o+9, o+10])
u.add_TopologyAttr('bonds', [tuple(x) for x in bonds])

lx = dx/2
lz = pbcz/2 - 65/2
coordinates = []
for i in range(Nx):
    for j in range(Nx):
        x = lx + dx*i
        y = lx + dx*j

        coordinates.append([x,y,lz])
        coordinates.append([x,y,lz+sigma])
        coordinates.append([x,y,lz+sigma*2])
        coordinates.append([x+0.5*sigma,y,lz+sigma*(2+0.866)])
        coordinates.append([x+0.5*sigma,y,lz+sigma*(3+0.866)])
        coordinates.append([x+0.5*sigma,y,lz+sigma*(4+0.866)])
        coordinates.append([x+0.5*sigma,y,lz+sigma*(5+0.866)])
        coordinates.append([x-0.5*sigma,y,lz+sigma*(2+0.866)])
        coordinates.append([x-0.5*sigma,y,lz+sigma*(3+0.866)])
        coordinates.append([x-0.5*sigma,y,lz+sigma*(4+0.866)])
        coordinates.append([x-0.5*sigma,y,lz+sigma*(5+0.866)])

ux = dx/2
uz = pbcz/2 + 65/2
for i in range(Nx):
    for j in range(Nx):
        x = ux + dx*i
        y = ux + dx*j

        coordinates.append([x,y,uz])
        coordinates.append([x,y,uz-sigma])
        coordinates.append([x,y,uz-sigma*2])
        coordinates.append([x+0.5*sigma,y,uz-sigma*(2+0.866)])
        coordinates.append([x+0.5*sigma,y,uz-sigma*(3+0.866)])
        coordinates.append([x+0.5*sigma,y,uz-sigma*(4+0.866)])
        coordinates.append([x+0.5*sigma,y,uz-sigma*(5+0.866)])
        coordinates.append([x-0.5*sigma,y,uz-sigma*(2+0.866)])
        coordinates.append([x-0.5*sigma,y,uz-sigma*(3+0.866)])
        coordinates.append([x-0.5*sigma,y,uz-sigma*(4+0.866)])
        coordinates.append([x-0.5*sigma,y,uz-sigma*(5+0.866)])

u.atoms.positions = coordinates
u.dimensions = [dx*(Nx), dx*(Nx), pbcz, 90, 90, 90]
u.atoms.write('membrane.gro')


with open('topol.top', 'w') as W:
    W.write('''#include "toppar/cg.itp"
#include "toppar/DOPC.itp"
#include "toppar/TRIO.itp"

[ system ]
; name
My own CG system

[ molecules ]
; name    n_residues
  DOPC    %d
''' %N)


