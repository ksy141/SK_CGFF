import numpy as np

### Set PGL - TGL coefficient
PT = 2.5


### Make a Gaussian table
rcut = 2.4
dr   = 0.3
delr = 0.002

r = np.arange(0, rcut + 1 + dr + delr/2, delr)
fr  = np.zeros(r.shape)
fr_ = np.zeros(r.shape)

gr = -np.exp(-2 * r**2)
gr_ = (-4 * r) * np.exp(-2 * r**2)

r[0] = 1
hr = 1/r**12
hr_ = 12 * 1/r**13
r[0] = 0
hr[r < 0.1] = 0
hr_[r < 0.1] = 0

np.savetxt('table.xvg', np.transpose([r, fr, fr_, gr, gr_, hr, hr_]), fmt='%8e')


### Make cg.itp
ofile = open('cg.itp', 'w')
ofile.write('''
[ defaults ]
; you should use 1 1 for user defined potentials
; nbfunc(1 for LJ 2 for Buckingham)   comb-rule
  1                                   1

[ atomtypes ]
; name mass charge ptype   C6     C12
  NC3  {NC3:6.1f}   0.0    A     0      0
  PO4  {PO4:6.1f}   0.0    A     0      0
  PGL  {PGL:6.1f}   0.0    A     0      0
  CT1  {CT1:6.1f}   0.0    A     0      0
  CT2  {CT2:6.1f}   0.0    A     0      0
  TGL  {TGL:6.1f}   0.0    A     0      0

[ nonbond_params ]
; i    j    func    C6     C12
'''.format(NC3=87, PO4=95, PGL=157, CT1=55.9, CT2=55.9, TGL=215))

######## Potential for nonbonded #########
# V(r) = k0 f(r) + C6 g(r) + C12 h(r)
# f(r) = 0
# g(r) = -exp(-2 r**2)
# h(r) = 1/r**12
# C6 = 0 for NC3, PO4, PGL
# C6 = 1.5 kcal/mol * g(r)
# C12 = 4 (e1e2)**0.5 * (0.5*(s1 + s2))**12
# e1 = e2 = 0.0028 kcal/mol

eps = 0.0028
kcal2kj = 4.184
A2nm = 0.1

from collections import namedtuple
Atom = namedtuple('Atom', 'name sigma epsilon C6')
NC3 = Atom('NC3', 9.5, eps, 0)
PO4 = Atom('PO4', 7.0, eps, 0)
PGL = Atom('PGL', 7.5, eps, 0)
CT1 = Atom('CT1', 6.8, eps, 1.5)
CT2 = Atom('CT2', 6.9, eps, 0.8)

TGL = Atom('TGL', 8.0, eps, 0)
Atoms = [NC3, PO4, PGL, CT1, CT2, TGL]
print(Atoms, '\n\n')

for i, ai in enumerate(Atoms):
    for j, aj in enumerate(Atoms):
        if i <= j:
            sij = 0.5 * (ai.sigma + aj.sigma) * A2nm
            eij = np.sqrt(ai.epsilon * aj.epsilon) * kcal2kj
            C12 = 4 * eij * sij**12

            if ai.name == 'CT1' and aj.name == 'CT2':
                C6 = 1.0 * kcal2kj

            elif ai.name == 'PGL' and aj.name == 'TGL':
                C6 = PT * kcal2kj

            elif ai.name == 'TGL' and aj.name == 'TGL':
                C6 = 1.5 * kcal2kj

            else:
                C6 = np.sqrt(ai.C6 * aj.C6) * kcal2kj

            ofile.write('{:5s} {:5s} {:<5d} {:8E}     {:8E}\n'.format(ai.name, aj.name, 1, C6, C12))
ofile.close()

