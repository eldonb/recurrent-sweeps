"""
sweep model
"""

import random
import numpy as np
import tskit
import msprime


Ne = 1e4
L = 10000000
# L = 30954429  # Length of simulated region
num_reps = 100
u = float(0)

# define hard sweep model
models = []
for i in range(100):
    u = random.uniform(0, 1) if i < 10 else random.uniform(0, 1)
    models.append(msprime.StandardCoalescent(duration=u))
    models.append(
        msprime.SweepGenicSelection(
        position=random.randint(1, L-1),  # middle of chrom
        start_frequency=1.0 / (2 * Ne),
        end_frequency=1.0 - (1.0 / (2 * Ne)),
        s=0.45,
        dt=1e-6,
        )
    )
models.append(msprime.StandardCoalescent())
print(models[:5])


reps = msprime.sim_ancestry(
    136,
    model=models,
    population_size=Ne,
    recombination_rate=1e-7,
    sequence_length=L,
    num_replicates=num_reps,
    random_seed=random.randint(15448, 55940303)
)

wins = np.linspace(0, L, 51)
mids = (wins[1:] + wins[:-1]) / 2

v = 0
n=136
afs=np.zeros(n-1, dtype=float)

ID = 'E'

for ts in reps:
    np.savetxt('branchespi_resout' + str(ID) + str(v),  ts.diversity(windows=wins, mode='branch'))
    mts=msprime.sim_mutations(ts, rate=1e-8,  random_seed=random.randint(1432, 4690933))
    m = mts.allele_frequency_spectrum([np.arange(n)],polarised=True, span_normalise=False)
    v = v + 1
    print(m.size)
    afs = afs + (m[1:136]/float(np.sum(m[1:136])))


print(afs/float(num_reps))
np.savetxt("recurrentsweeps_sfs_resout" + str(ID), afs/float(num_reps))
