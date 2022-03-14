#!~/.conda/envs/lra/bin/python3
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import re
import itertools
from ortools.sat.python import cp_model
import Levenshtein as lv
import multiprocessing 
from statistics import mean

os.chdir(r"/project/mchaisso_100/cmb-16/jingwenr/trfCall")

if len(sys.argv) == 2:
    input_file = sys.argv[1]
else:
    print("require file_containing_motifs_and_counts")
    exit(1)

class MotifsInfo:
    def __init__(self):
        self.motifs = []
        self.motifs_counts = []
        self.editdist = []
        self.numOfVNTR = 0
        self.delta = []
        self.coords = []
      
    """
    read the 2nd and 3rd column of the file (ranked_motifs, ranked_motif_counts)
    each row is a vntr site
    """
    def readFile (self, input_file):
        with open(input_file, 'r') as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                if idx == 0: #### Update the idx constraint
                    continue
                fields = re.split('\t', line.strip('\n'))
                self.coords.append(fields[0])
                self.motifs.append(fields[1].split(','))
                self.motifs_counts.append([int(cnt) for cnt in fields[2].split(',')])
        self.numOfVNTR = len(self.motifs)
        assert (len(self.motifs) == len(self.motifs_counts)),"motifs and motifs_counts have different lengths!"
                
    def getMotifsForOneVNTR(self, i):
        return self.motifs[i]

    def getMotifsCountsForOneVNTR(self, i):
        return self.motifs_counts[i]
        
    """
    compute the pairwise Levenshtein edit distance between any pair of motifs
    """ 
    def pairwiseLevenshteinDist(self):
        for i in range(self.numOfVNTR):
            numOfMotifs = len(self.motifs[i])
            dist = [[0 for h in range(numOfMotifs)] for k in range(numOfMotifs)]
            for k in range(numOfMotifs):
                for h in range(numOfMotifs):
                    if k <= h: 
                        continue
                    dist[k][h] = lv.distance(self.motifs[i][k], self.motifs[i][h])
                    dist[h][k] = dist[k][h]
            self.editdist.append(dist)
    
    """
    compute the delta array for each VNTR locus [numOfmotifs * 0.2 quantile edit distance, ..., numOfmotifs * 0.8 quantile edit distance]
    """
    def decideDelta(self):
        self.delta = [[] for i in range(self.numOfVNTR)]
        for idx_vntr in range(self.numOfVNTR):
            numOfmotifs = len(self.editdist[idx_vntr])
            temp = np.array(list(itertools.chain(*self.editdist[idx_vntr])))
            self.delta[idx_vntr] = [int(np.quantile(temp, i / 10) * numOfmotifs) for i in range(1, 10)]
            print(self.delta[idx_vntr])
            # l = int(np.quantile(temp, 0.1))
            # u = int(np.quantile(temp, 0.9))
            # intv = ((u - l) *  numOfmotifs) // 9
            # self.delta[idx_vntr] = range(l * numOfmotifs, u * numOfmotifs + 1, intv)
    """
    plot out the histogram of edit distance for a vntr site
    """
    def HisteditDist(self):
        for idx_vntr in range(self.numOfVNTR):
            dist = np.array(list(itertools.chain(*self.editdist[idx_vntr])))
            fig = plt.hist(dist, density=False, bins=100)  # density=False would make counts
            plt.ylabel('Count')
            plt.xlabel('Edit Distance')
            plt.title("Histogram of Pairwise Edit Distance")
            plt.show()       
            plt.savefig("hist/vntr_motif_editdist." + str(idx_vntr) + ".png")
            plt.close()

"""
build up the ILP solver
"""
def ILPsolver(motifs, motifs_counts, costs, idx_vntr, delta):
    model = cp_model.CpModel()
    num_originalMotifs, num_afterMotifs = len(costs), len(costs)

    """
    create the variables
    """
    x = []
    for i in range(num_originalMotifs):
        t = []
        for j in range(num_afterMotifs):
            t.append(model.NewBoolVar(f'x[{i},{j}]'))
        x.append(t)

    y = []
    for j in range(num_afterMotifs):
        y.append(model.NewBoolVar(f'y[{j}]'))

    """
    create the constraints
    """
    # Each motif is mapped to one motif
    for i in range(num_originalMotifs):
        model.Add(sum(x[i][j] for j in range(num_afterMotifs)) == 1)

    # The upper bound
    for j in range(num_afterMotifs):
        model.Add(1 - sum(x[i][j] for i in range(num_originalMotifs)) <= 1 - y[j])

    # The lower bound
    for j in range(num_afterMotifs):
        model.Add(1 - sum(x[i][j] for i in range(num_originalMotifs)) >= -num_afterMotifs * y[j] + 1)

    # The total cost is less than or equal to \delta
    model.Add(sum(motifs_counts[i] * sum(x[i][j] * int(costs[i][j]) for j in range(num_afterMotifs)) for i in range(num_originalMotifs)) <= delta)

    """
    create objective function
    """
    objective_terms = []
    for j in range(num_afterMotifs):
        objective_terms.append(y[j])
    model.Minimize(sum(objective_terms))

    """
    invoke the solver
    """
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    """
    print out the solution
    """
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print()
        print(f'Delta = {delta}')
        print(f'Original Domain Size = {num_originalMotifs}')
        print(f'Domain size = {int(solver.ObjectiveValue())}')
        print()

        # for i in range(num_originalMotifs):
        #     for j in range(num_afterMotifs):
        #         if solver.BooleanValue(x[i][j]):
        #             print(f'motif {i} assigned to motif {j} Cost = {costs[i][j]}')
    else:
        print('No solution found.')

    return num_originalMotifs, int(solver.ObjectiveValue())

"""
plot out original_motifs_size VS after_motifs_size over the delta array for each vntr
"""
def PlotMotifSizeChangeWithDelta(original_domain_size, after_domain_size, delta, idx_vntr):
    plt.plot(delta, original_domain_size, label = 'original motifs size', color = 'blue', marker = 'o')
    plt.plot(delta, after_domain_size, label = 'efficient motifs size', color = 'red', marker = '^')
    plt.legend()
    plt.ylabel('motif size')
    plt.xlabel('delta')
    plt.title("efficient motif size change with delta")
    plt.show()
    plt.savefig("png/EfficientMotifSizeChangeWithDelta." + str(idx_vntr) + ".png", dpi=300)
    plt.close()    

"""
multiple processing
"""
motifsInfo = MotifsInfo() 
motifsInfo.readFile(input_file)
motifsInfo.pairwiseLevenshteinDist()
motifsInfo.decideDelta()
original_motifs_sizes = [[] for i in range(motifsInfo.numOfVNTR)]
after_motifs_sizes = [[] for i in range(motifsInfo.numOfVNTR)]
print("finish reading files")
    
def mp_worker(idx_vntr):

    original_domain_size = []
    after_domain_size = []
    for delta in motifsInfo.delta[idx_vntr]:
        print(delta)
        original, after = ILPsolver(motifsInfo.motifs[idx_vntr], motifsInfo.motifs_counts[idx_vntr], motifsInfo.editdist[idx_vntr], idx_vntr, delta)
        original_domain_size.append(original)
        after_domain_size.append(after)
    return idx_vntr, original_domain_size, after_domain_size

def mp_handler():
    p = multiprocessing.Pool(2)
    all_vntr = range(motifsInfo.numOfVNTR)
    for idx_vntr, original_domain_size, after_domain_size in p.imap(mp_worker, all_vntr):
        print("outputing " + str(idx_vntr))
        original_motifs_sizes[idx_vntr] = original_domain_size
        after_motifs_sizes[idx_vntr] = after_domain_size
        PlotMotifSizeChangeWithDelta(original_domain_size, after_domain_size, motifsInfo.delta[idx_vntr], idx_vntr)

if __name__=='__main__':
    mp_handler()

