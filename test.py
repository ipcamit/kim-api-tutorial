from ase.build import bulk
from ase.calculators.kim.kim import KIM
si = bulk("Si")
model = KIM("LJSi_MO_111111111110_000")
si.calc = model
print(si.get_potential_energy())
