import MDAnalysis as mda
import numpy as np

### N = 3600 * 2; dN = 180 * 2
### N = 1600 * 2; dN = 80 * 2
pbcz  = 200
N     = 900 * 2
dN    = 45 * 2
Nx    = int(np.sqrt(N/2))
dx    = 8

sigma  = 5 #A
PL_natoms = 11
TG_natoms = 13

n_atoms = (N-dN)*PL_natoms + dN*TG_natoms
n_residues = N
PL_resindex = np.repeat(np.arange(0, N-dN), PL_natoms)
TG_resindex = np.repeat(np.arange(N-dN, N), TG_natoms)
atom_resindex = np.append(PL_resindex, TG_resindex)

u = mda.Universe.empty(n_atoms = n_atoms,
                       n_residues = n_residues,
                       atom_resindex = atom_resindex,
                       trajectory = True)

u.add_TopologyAttr('name', ['NC3','PO4','PGL',
                            'C21', 'C22', 'C23', 'C24',
                            'C31', 'C32', 'C33', 'C34'] * (N-dN) +
                           ['TGL',
                            'C11', 'C12', 'C13', 'C14',
                            'C21', 'C22', 'C23', 'C24',
                            'C31', 'C32', 'C33', 'C34'] * dN
                            )

u.add_TopologyAttr('type', ['NC3','PO4','PGL',
                            'CT1','CT1','CT2','CT2',
                            'CT1','CT1','CT2','CT2'] * (N-dN) +
                            ['TGL',
                             'CT1', 'CT2', 'CT2', 'CT2',
                             'CT1', 'CT2', 'CT2', 'CT2',
                             'CT1', 'CT2', 'CT2', 'CT2'] * dN)

u.add_TopologyAttr('resname', ['DOPC'] * (N-dN) + ['TRIO'] * dN)
u.add_TopologyAttr('resid', np.arange(1, N+1))

bonds = []
for o in range(0, (N-dN)*11, 11):
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

for o in range((N-dN)*11, (N-dN)*11 + dN*13, 13):
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

TG_locs = np.arange(N/2, dtype=np.int)
np.random.shuffle(TG_locs)

TG_up_coords = []
lx = dx/2
lz = pbcz/2 - 65/2
coordinates = []
for i in range(Nx):
    for j in range(Nx):
        x = lx + dx*i
        y = lx + dx*j

        if (i*Nx + j) in TG_locs[:int(dN/3)]:
            TG_up_coords.append([x, y])
            continue

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


TG_locs = np.arange(N/2, dtype=np.int)
np.random.shuffle(TG_locs)
TG_dw_coords = []
ux = dx/2
uz = pbcz/2 + 65/2
for i in range(Nx):
    for j in range(Nx):
        x = ux + dx*i
        y = ux + dx*j

        if (i*Nx + j) in TG_locs[:int(2*dN/3)]:
            TG_dw_coords.append([x, y])
            continue

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

for TG_coord in TG_up_coords:
    x, y = TG_coord
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

for TG_coord in TG_dw_coords:
    x, y = TG_coord
    coordinates.append([x,y,uz-sigma*2])
    coordinates.append([x-0.5*sigma,y-0.2887*sigma,uz-sigma*(2+0.8165)])
    coordinates.append([x-0.5*sigma,y-0.2887*sigma,uz-sigma*(3+0.8165)])
    coordinates.append([x-0.5*sigma,y-0.2887*sigma,uz-sigma*(4+0.8165)])
    coordinates.append([x-0.5*sigma,y-0.2887*sigma,uz-sigma*(5+0.8165)])
    coordinates.append([x+0.5*sigma,y-0.2887*sigma,uz-sigma*(2+0.8165)])
    coordinates.append([x+0.5*sigma,y-0.2887*sigma,uz-sigma*(3+0.8165)])
    coordinates.append([x+0.5*sigma,y-0.2887*sigma,uz-sigma*(4+0.8165)])
    coordinates.append([x+0.5*sigma,y-0.2887*sigma,uz-sigma*(5+0.8165)])
    coordinates.append([x+0.0*sigma,y+0.5774*sigma,uz-sigma*(2+0.8165)])
    coordinates.append([x+0.0*sigma,y+0.5774*sigma,uz-sigma*(3+0.8165)])
    coordinates.append([x+0.0*sigma,y+0.5774*sigma,uz-sigma*(4+0.8165)])
    coordinates.append([x+0.0*sigma,y+0.5774*sigma,uz-sigma*(5+0.8165)])

u.atoms.positions = coordinates
u.dimensions = [dx*(Nx), dx*(Nx), pbcz, 90, 90, 90]
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
''' %(N-dN, dN))
