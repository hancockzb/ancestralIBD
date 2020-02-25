cd ~/Desktop
source activate msprime-env
python3
import msprime, pyslim
import allel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
np.set_printoptions(threshold=np.inf)
np.set_printoptions(formatter={'float': lambda x: "{0:0.6f}".format(x)})

#Start a population with a coalescent history
ts = msprime.simulate(sample_size=2000, Ne=4000, mutation_rate=0.0, recombination_rate=1e-8)
ts = pyslim.annotate_defaults(ts, model_type="WF", slim_generation=1)
ts.dump("msprime_coal.trees")

#Move to SLiM:

./slim mig1_50_slim.txt
./slim mig01_50_slim.txt
./slim mig001_50_slim.txt
./slim mig1_25_slim.txt
./slim mig01_25_slim.txt
./slim mig001_25_slim.txt
./slim mig1_10_slim.txt
./slim mig01_10_slim.txt
./slim mig001_10_slim.txt
./slim mig1_5_slim.txt
./slim mig01_5_slim.txt
./slim mig001_5_slim.txt
./slim mig1_1_slim.txt
./slim mig01_1_slim.txt
./slim mig001_1_slim.txt

#Overlay neutral mutations

ts = pyslim.load("demo_100_1.trees").simplify()
mutated = msprime.mutate(ts, rate=1e-7, model=msprime.InfiniteSites(msprime.NUCLEOTIDES), random_seed=1)
mutated.dump("mig1_10_mut.trees")
ts = pyslim.load("mig1_50_mut.trees").simplify()
ts.num_sites

####Methods for paper!####

#Pairwise divergence across the genome, take a sample from each subpop, write to file
A = ts.samples()[2000:3999]
B = ts.samples()[4000:5999]
windows = np.linspace(0, ts.sequence_length, num=1000)
#or
windows = ts.breakpoints()
#Sequence divergence
x = ts.divergence([A,B], windows=windows)
a = str(x)
f = open("div_45_mig001_1.txt", "w")
f.write(a)
f.close()

#Get TMRCA across genome for subsample and write to file
samples09 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
             6000, 6005, 6010, 6015, 6020, 6025, 6030, 6035, 6040, 6045]
samples45 = [2000, 2005, 2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045,
           4000, 4005, 4010, 4015, 4020, 4025, 4030, 4035, 4040, 4045]
ts_sub = ts.simplify(samples45)
tmrca = np.zeros(ts.num_trees)
breakpoints = np.zeros(ts.num_trees)
for tree in ts.trees():
    tmrca[tree.index] = tree.time(tree.root)
    breakpoints[tree.index] = tree.interval[0]
a = str([tmrca, breakpoints])
f = open("demo_100_1.txt", "w")
f.write(a)
f.close()
plt.ylabel("T_mrca (Generations)") 
plt.xlabel("Position (kb)") 
plt.plot(breakpoints / 1000, tmrca, "o");

#IBD tract length

def ibd_segments(ts, a, b): 
	trees_iter = ts.trees() 
	tree = next(trees_iter) 
	last_mrca = tree.mrca(a, b) 
	last_left = 0 
	segment_lengths = []
	for tree in trees_iter: 
		mrca = tree.mrca(a, b) 
		if mrca != last_mrca:
			left = tree.interval[0] 
			segment_lengths.append(left - last_left) 
			last_mrca = mrca
			last_left = left
	segment_lengths.append(ts.sequence_length - last_left) 
	return np.array(segment_lengths) / ts.sequence_length
print(ibd_segments(ts, 0, 1))
print(ibd_segments(ts, 0, 6000))
sns.distplot(ibd_segments(ts, a, b), label="Within population")
sns.distplot(ibd_segments(ts, a, b), label="Between populations") 
plt.xlim(-0.0001, 0.02)
plt.legend()
plt.xlabel("Fraction of genome length")
plt.ylabel("Count")
plt.show()

#Subsample individuals, output variant haplotypes
samples09 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
             6000, 6005, 6010, 6015, 6020, 6025, 6030, 6035, 6040, 6045,
             8000]
ts_010 = ts.simplify(samples09)
samples45 = [2000, 2005, 2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045,
           4000, 4005, 4010, 4015, 4020, 4025, 4030, 4035, 4040, 4045,
           8000]
ts_45 = ts.simplify(samples45)
haps = []
for i in ts_45.haplotypes():
    haps.append(i)
sequence_IDs = []
for i in range(len(haps)):
    sequence_IDs.append(f'sample_{ts_45.samples()[i]}_pop_{ts_45.node(i).population}')
with open('ts_45_mig01_5.fas', 'w') as f:
    for i in range(len(haps)):
        f.write(f'>{sequence_IDs[i]}\n{haps[i]}\n')
