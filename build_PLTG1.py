import MDAnalysis as mda
import numpy as np

### N = 3600 * 2; dN = 180 * 2
### N = 1600 * 2; dN = 80 * 2
pbcz  = 100
N     = 64 * 2
dN    = 1
Nx    = int(np.sqrt(N/2))
dx    = 8

sigma  = 5 #A
PL_natoms = 11
TG_natoms = 13

n_atoms = N*PL_natoms + dN*TG_natoms
n_residues = N + dN
PL_resindex = np.repeat(np.arange(0, N), PL_natoms)
TG_resindex = np.repeat(np.arange(N, N+dN), TG_natoms)
atom_resindex = np.append(PL_resindex, TG_resindex)

u = mda.Universe.empty(n_atoms = n_atoms,
                       n_residues = n_residues,
                       atom_resindex = atom_resindex,
                       trajectory = True)

u.add_TopologyAttr('name', ['NC3','PO4','PGL',
                            'C21', 'C22', 'C23', 'C24',
                            'C31', 'C32', 'C33', 'C34'] * N +
                           ['TGL',
                            'C11', 'C12', 'C13', 'C14',
                            'C21', 'C22', 'C23', 'C24',
                            'C31', 'C32', 'C33', 'C34'] * dN
                            )

u.add_TopologyAttr('type', ['NC3','PO4','PGL',
                            'CT1','CT1','CT2','CT2',
                            'CT1','CT1','CT2','CT2'] * N +
                            ['TGL',
                             'CT1', 'CT2', 'CT2', 'CT2',
                             'CT1', 'CT2', 'CT2', 'CT2',
                             'CT1', 'CT2', 'CT2', 'CT2'] * dN)

u.add_TopologyAttr('resname', ['DOPC'] * N + ['TRIO'] * dN)
u.add_TopologyAttr('resid', np.arange(1, N+dN+1))

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

for o in range(N*11, N*11 + dN*13, 13):
    bonds.append([o+0, o+1])
    bonds.append([o+0, o+5])
    bonds.append([o+0, o+9])
    bonds.append([o+1, o+2])
    bonds.append([o+2, o+3])
    bonds.append([o+3, o+4])
    bonds.append([o+5, o+6])
    bonds.append([o+6, o+7])
    bonds.append([o+7, o+8])
    bonds.append([o+9, o+10])
    bonds.append([o+10, o+11])
    bonds.append([o+11, o+12])

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


x = ux + dx * Nx
y = ux + dx * Nx
coordinates.append([x,y,lz+sigma*2])
coordinates.append([x-0.5*sigma,y-0.2887*sigma,lz+sigma*(2+0.8165)])
coordinates.append([x-0.5*sigma,y-0.2887*sigma,lz+sigma*(3+0.8165)])
coordinates.append([x-0.5*sigma,y-0.2887*sigma,lz+sigma*(4+0.8165)])
coordinates.append([x-0.5*sigma,y-0.2887*sigma,lz+sigma*(5+0.8165)])
coordinates.append([x+0.5*sigma,y-0.2887*sigma,lz+sigma*(2+0.8165)])
coordinates.append([x+0.5*sigma,y-0.2887*sigma,lz+sigma*(3+0.8165)])
coordinates.append([x+0.5*sigma,y-0.2887*sigma,lz+sigma*(4+0.8165)])
coordinates.append([x+0.5*sigma,y-0.2887*sigma,lz+sigma*(5+0.8165)])
coordinates.append([x+0.0*sigma,y+0.5774*sigma,lz+sigma*(2+0.8165)])
coordinates.append([x+0.0*sigma,y+0.5774*sigma,lz+sigma*(3+0.8165)])
coordinates.append([x+0.0*sigma,y+0.5774*sigma,lz+sigma*(4+0.8165)])
coordinates.append([x+0.0*sigma,y+0.5774*sigma,lz+sigma*(5+0.8165)])

u.atoms.positions = coordinates
u.dimensions = [dx*(Nx+1), dx*(Nx+1), pbcz, 90, 90, 90]
u.atoms.write('step5.gro')

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
  TRIO    %d
''' %(N, dN))



