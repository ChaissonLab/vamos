#!~/.conda/envs/lra/bin/python3
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import re
import itertools
from ortools.sat.python import cp_model
import Levenshtein as lv
# import multiprocessing 
from ray.util.multiprocessing import Pool
from statistics import mean
import random 
import traceback
import time 

#os.chdir(r"/project/mchaisso_100/cmb-16/jingwenr/trfCall")
random.seed(0)

if len(sys.argv) == 5:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    numOfProcessors = int(sys.argv[3])
    delta_threshold = float(sys.argv[4])
else:
    print("require file_containing_motifs_and_counts output_file numOfProcessors delta_threshold")
    exit(1)

print(delta_threshold)

class MotifsInfo:
    def __init__(self):
        self.motifs = []
        self.motifs_counts = []
        self.editdist = []
        self.numOfVNTR = 0
        self.delta = []
        self.coords = []
        self.total_motifs_cnt = []
      
    """
    read the 2nd and 3rd column of the file (ranked_motifs, ranked_motif_counts)
    each row is a vntr site
    """
    def readFile (self, input_file):
        with open(input_file, 'r') as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                # if idx == 0: #### Update the idx constraint
                #     continue
                fields = re.split('\t', line.strip('\n'))
                # if fields[0] == "chr6_157310355-157314362" or fields[0] == "chr7_158146914-158149196": continue
                self.coords.append(fields[0])
                self.motifs.append(fields[1].split(','))
                self.motifs_counts.append([int(cnt) for cnt in fields[2].split(',')])
                self.total_motifs_cnt.append(sum(self.motifs_counts[-1])) ### new added
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
    compute the delta based on the mean pairwise Levenshtein edit distance
    """
    def decideDelta(self):
        self.delta = [0] * self.numOfVNTR
        for idx_vntr in range(self.numOfVNTR):
            numOfmotifs = len(self.editdist[idx_vntr])
            if numOfmotifs == 1:
                self.delta[idx_vntr] = 0 
                continue

            temp = []
            for k in range(numOfmotifs):
                for h in range(numOfmotifs):
                    if k < h: # k replace h                        
                        temp.extend([self.editdist[idx_vntr][k][h]] * self.motifs_counts[idx_vntr][k])
                    elif h < k:
                        temp.extend([self.editdist[idx_vntr][h][k]] * self.motifs_counts[idx_vntr][h])
            temp_arr = np.array(temp)
            temp_arr.sort()
            # self.delta[idx_vntr] = int(np.quantile(temp_arr, delta_threshold) * self.total_motifs_cnt[idx_vntr]) # delta strategy 1
            # self.delta[idx_vntr] = int(np.quantile(temp_arr, delta_threshold) * int(self.total_motifs_cnt[idx_vntr] * delta_threshold)) # delta strategy 2
            self.delta[idx_vntr] = int(np.quantile(temp_arr, delta_threshold) * int(self.total_motifs_cnt[idx_vntr] * delta_threshold)) # delta strategy 3
        return 


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
    start = time.time()
    print("processing VNTR " + str(idx_vntr))
    model = cp_model.CpModel()
    num_originalMotifs, num_afterMotifs = len(costs), len(costs)
    assert(len(motifs) == len(costs)), "motifs and costs have different lengths!"

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

    # Declare our intermediate boolean variable.
    b = []
    for i in range(num_originalMotifs):
        t = []
        for j in range(num_afterMotifs):
            t.append(model.NewBoolVar(f'b[{i},{j}]'))
        b.append(t)


    # Implement b[i][j] == (motifs_counts[i] < motifs_counts[j]).
    # for i in range(num_originalMotifs):
    #     for j in range(num_afterMotifs):
    #         model.Add(motifs_counts[i] < motifs_counts[j]).OnlyEnforceIf(b[i][j])
    #         model.Add(motifs_counts[i] >= motifs_counts[j]).OnlyEnforceIf(b[i][j].Not())

    for i in range(num_originalMotifs):
        for j in range(num_afterMotifs):
            model.Add(motifs_counts[i] > motifs_counts[j]).OnlyEnforceIf(b[i][j])
            model.Add(motifs_counts[i] <= motifs_counts[j]).OnlyEnforceIf(b[i][j].Not())


    for i in range(num_originalMotifs):
        for j in range(num_afterMotifs):
            model.Add(x[i][j] == 0).OnlyEnforceIf(b[i][j])
            model.Add(x[i][j] >= 0).OnlyEnforceIf(b[i][j].Not())



    # # The total cost is less than or equal to \delta
    model.Add(sum(motifs_counts[i] * sum(x[i][j] * int(costs[i][j]) for j in range(num_afterMotifs)) for i in range(num_originalMotifs)) <= delta)

    # """
    # create objective function
    # """
    # objective_terms = []
    # for j in range(num_afterMotifs):
    #     objective_terms.append(y[j])
    # model.Minimize(sum(objective_terms))
    
    # The total efficent motifs <= 255
    # model.Add(sum(y) <= 255)

    """
    create objective function
    minimizing the sum of efficient motif set size and replacement cost
    """
    # model.Minimize(sum(y))

    model.Minimize(delta * sum(y) + sum(motifs_counts[i] * sum(x[i][j] * int(costs[i][j]) for j in range(num_afterMotifs)) for i in range(num_originalMotifs)))

    """
    invoke the solver
    """
    solver = cp_model.CpSolver()
    print("invoke the solver")

     # Sets a time limit of 300 seconds
    solver.parameters.max_time_in_seconds = 600.0
    solver.parameters.num_search_workers = 2
    # solver.parameters.enumerate_all_solutions = True

    print("start the solver")
    status = solver.Solve(model)
    print("check the status")

    """
    print out the solution
    """
    aftermotifs = []
    aftercounts = [0] * num_originalMotifs
    mapcosts = []
    clean_motifs = []
    clean_motifs_counts = []
    status_str = ''
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print()
        if status == cp_model.OPTIMAL:
            status_str = "OPTIMAL"
            print(f'OPTIMAL solution found')
        else:
            status_str = "FEASIBLE"
            print(f'A solution found, but may not be optimal')

        print(f'Original Domain Size = {num_originalMotifs}')
        print(f'Domain size = {int(solver.ObjectiveValue())}')
        print()

        for i in range(num_originalMotifs):
            check = False
            for j in range(num_afterMotifs):
                if solver.BooleanValue(x[i][j]):
                    aftermotifs.append(motifs[j])
                    aftercounts[j] += motifs_counts[i]
                    mapcosts.append(costs[i][j])
                    check = True
                    # print(f'motif {i} assigned to motif {j} Cost = {costs[i][j]}')
                    break
            assert check == True, f'motif i - {i} is mapped to nothing!'
        assert len(motifs) == len(aftermotifs), f'motifs - {len(motifs)} and aftermotifs - {len(aftermotifs)} have different lengths!'
    else:
        status_str = "INFEASIBLE"
        print('No solution found.')
        return idx_vntr, list(zip(motifs, motifs)), ['0' for cnt in range(num_originalMotifs)], num_originalMotifs, num_originalMotifs, \
            motifs, motifs_counts, ['1' for cnt in motifs_counts], delta, motifs_counts, status_str
        # return idx_vntr, [], mapcosts, num_originalMotifs, int(solver.ObjectiveValue()), clean_motifs, clean_motifs_counts, [], delta, motifs_counts

    # Statistics.
    print('\nStatistics')
    print(f'  VNTR idx : {idx_vntr}')
    print(f'  status   : {solver.StatusName(status)}')
    print(f'  conflicts: {solver.NumConflicts()}')
    print(f'  branches : {solver.NumBranches()}')
    print(f'  wall time: {solver.WallTime()} s')
    
    for i in range(len(aftercounts)):
    	if aftercounts[i] > 0:
    		clean_motifs.append(motifs[i])
    		clean_motifs_counts.append(aftercounts[i])

    assert len(motifs) == len(aftermotifs), f'motifs - {len(motifs)} and aftermotifs - {len(aftermotifs)} have different lengths!'

    indicator = ['0' if cnt == 0 else '1' for cnt in aftercounts]
    end = time.time()
    print(f'finish vntr-{idx_vntr} (time: {(end - start)/60} min)')
    
    # return idx_vntr, list(zip(motifs, aftermotifs)), aftercounts, mapcosts, num_originalMotifs, int(solver.ObjectiveValue())
    return idx_vntr, list(zip(motifs, aftermotifs)), mapcosts, num_originalMotifs, int(solver.ObjectiveValue()), \
        clean_motifs, clean_motifs_counts, indicator, delta, motifs_counts, status_str

"""
plot out the histogram of mapping cost for each vntr site
"""
def Histmapcost(idx_vntr, mapcosts, _mean):
    fig = plt.hist(mapcosts, density=False)  # density=False would make counts
    plt.axvline(_mean, color='red', linestyle='dashed', linewidth=2)    
    plt.ylabel('Count')
    plt.xlabel('mapping cost')
    plt.title("Histogram of mapping cost")
    plt.show()       
    plt.savefig("mapcost/vntr_motif_mapcost." + str(idx_vntr) + ".png", dpi=300)
    plt.close()

"""
plot out the histogram of mapped motif size for all the VNTR loci
"""
def HistMappedMotifSize(after_domain_size):
    fig = plt.hist(after_domain_size, density=False, bin=100)  # density=False would make counts
    plt.ylabel('Count')
    plt.xlabel('Mapped/Efficient Motif Size')
    plt.title("Histogram off Mapped/Efficient Motif Size")
    plt.show()       
    plt.savefig("hist_mapped_motif_size.png", dpi=300)
    plt.close()

"""
plot out the original_motifs_size VS after_motifs_size
"""
def PlotMotifSizeChange(original_motifs_size, after_motifs_size):
	x = range(1, len(original_motifs_size) + 1)
	plt.plot(x, original_motifs_size, label = 'original motifs size', color = 'blue')
	plt.plot(x, after_motifs_size, label = 'mapped motifs size', color = 'red')
	plt.ylabel('motif size')
	plt.xlabel('vntr site')
	plt.legend()
	plt.show()
	plt.savefig("motif_size_comparison.png", dpi=300)
	plt.close()

"""
multiple processing
"""
start = time.time()
motifsInfo = MotifsInfo() 
motifsInfo.readFile(input_file)
motifsInfo.pairwiseLevenshteinDist()
motifsInfo.decideDelta()
original_motifs_size = []
after_motifs_size = []
end = time.time()
# motifsInfo.HisteditDist()
print(f'finish reading files (time: {(end - start) / 60} min)')
    
def mp_worker(idx_vntr):
    print("start ILP solver")
    try:
        return ILPsolver(motifsInfo.motifs[idx_vntr], motifsInfo.motifs_counts[idx_vntr], motifsInfo.editdist[idx_vntr], idx_vntr, motifsInfo.delta[idx_vntr])
    except Exception as e:
        print('Caught exception in worker thread')  
        traceback.print_exc() 
        print()
        raise e

def mp_handler():
    print("initiate multiprocessing")
    # pool = multiprocessing.Pool(numOfProcessors)
    pool = Pool(numOfProcessors)
    all_vntr = range(motifsInfo.numOfVNTR)
    print(all_vntr)
    print("start multiprocessing")

    with open(output_file, 'w') as f:
        # f.write('coordinate' + '\t'	+ 'orginal_motifs' + '\t' + 'mapped_motifs' + '\t' + 
        # 	    'mapping_cost' + '\t' + 'num_original_motifs' + '\t' + 'num_mapped_motifs' + '\t' + 
        # 	    'clean_motifs' + '\t' + 'clean_motifs_counts' + '\t' + 'indicator' + '\n')

        try:
            for idx_vntr, motifs_pair, mapcosts, original_domain_size, after_domain_size, clean_motifs, clean_motifs_counts, indicator, delta, motifs_counts, status in pool.imap_unordered(mp_worker, all_vntr):
                print("summarizing VNTR " + str(idx_vntr))
                if len(motifs_pair) > 0:
                    # Histmapcost(idx_vntr, mapcosts, motifsInfo.delta[idx_vntr] / len(motifsInfo.motifs[idx_vntr]))
                    original_motifs_size.append(original_domain_size)
                    after_motifs_size.append(after_domain_size)
                    f.write(motifsInfo.coords[idx_vntr] + '\t' + 
                            ','.join([str(motifs_pair[i][0]) for i in range(len(motifs_pair))]) + '\t' + 
                            ','.join([str(motifs_pair[i][1]) for i in range(len(motifs_pair))]) + '\t' + 
                            ','.join([str(cost) for cost in mapcosts]) + '\t' + 
                            str(original_domain_size) + '\t' + 
                            str(after_domain_size) + '\t' + 
                            ','.join(clean_motifs) + '\t' + 
                            ','.join([str(cnt) for cnt in clean_motifs_counts]) + '\t' + 
                            ','.join(indicator) + '\t' +
                            str(delta) + '\t' +
                            ','.join([str(cnt) for cnt in motifs_counts]) + '\t' + 
                            status + '\n')
        except:
            print ("Outer exception caught!")
    pool.close()
    pool.join()

if __name__=='__main__':
    mp_handler()
    # PlotMotifSizeChange(original_motifs_size, after_motifs_size)
    # HistMappedMotifSize(after_motifs_size)
