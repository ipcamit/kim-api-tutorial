from ase.build import bulk
from ase.calculators.kim.kim import KIM
si = bulk("Si")
model = KIM("LJSi_MO_111111111110_000")
si.calc = model
energy = si.get_potential_energy()
print(f"Energy: {energy:.6f} eV")
print(si.get_forces())
