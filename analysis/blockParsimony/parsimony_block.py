import sys
import datetime
import logging
import dendropy
import copy

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)


inBlocks = sys.argv[1]
inMappedNewick = sys.argv[2]
outTrees = sys.argv[3]

# read all selected block loci
logging.info(f'reading input block loci...')
blockTreeMap, treeDict, blockDict = {}, {}, {}
with open(inBlocks) as f:
    for line in f:
        fields = line.strip().split()
        chr,start,end,startT,endT,period,type,annotated,spanned,cons,blocks = fields[:11]
        alleles = [ int(a) for a in fields[11:] ]
        blockTreeMap[(chr,start,end)] = [startT,endT]
        treeDict[(chr,startT,endT)] = ''
        blockDict[(chr,start,end,cons)] = alleles

# read all trees
logging.info(f'reading input trees...')
with open(inMappedNewick) as f:
    for line in f:
        chr,_,_,_,startT,endT = line.strip().split()[:6]
        if (chr,startT,endT) in treeDict:
            treeDict[(chr,startT,endT)] = line.strip().split()[-1]
            #treeDict[(chr,startT,endT)] = dendropy.Tree.get(data=line.strip().split()[-1],schema='newick')


logging.info(f'inferring block history...')
out = open(outTrees, 'w')
for (chr,start,end,cons),alleles in blockDict.items():
    logging.info(f'locus: {chr}:{start}-{end};{cons}')
    startT,endT = blockTreeMap[(chr,start,end)]
    alleles = blockDict[(chr,start,end,cons)]
    print(treeDict[(chr,startT,endT)])
    tree = dendropy.Tree.get(data=treeDict[(chr,startT,endT)],schema='newick')
    # tree configuring
    #print(alleles)
    #print(tree.as_string("newick")) #tree.print_plot()
    #node = tree.find_node_with_taxon_label('0'); print(node)
    for i,node in enumerate(tree.postorder_node_iter()): node.label = i
    numNodeTotal = i+1
    logging.info(f'total number of nodes in the tree: {numNodeTotal}')
    #for i,node in enumerate(tree.postorder_node_iter()): print(node)


    ####################################################
    ##### dynamic programming to get the best tree #####
    ####################################################

    #### distance vector and traceback vector for all node
    # find the shortest and longest allele and the total number of possible states
    lenMin, lenMax = min(alleles), max(alleles)
    numStateTotal = lenMax-lenMin+1
    allStates = list(range(lenMin,lenMax+1,1))
    logging.info(f'total number of possible states: {numStateTotal}')
    logging.info(f'all possible states: {",".join(map(str,allStates))}')
    # the distance vector records the best distance of each state of each node
    distances = [ [float('inf')]*numStateTotal for i in tree.postorder_node_iter() ]
    # the traceback vector records the best state of the two child node for each current node state of a internal node
    traceback = [ [[0,0]]*numStateTotal for i in tree.postorder_node_iter() ]


    # assign initial scores for each leaf node by the allele feature
    for i,node in enumerate(tree.postorder_node_iter()):
        if node.is_leaf():
            sampleIdx = str(node.taxon).replace("'","")
            length = alleles[int(sampleIdx)]
            distances[i][length-lenMin] = 0
        else:
            pass

    # propagate the distance matrix by postorder (child comes before parent) for each internal node
    for i,node in enumerate(tree.postorder_node_iter()):
        if node.is_leaf():
            pass
            #logging.info(f'passing leaf node {i}.')
        else:
            child1, child2 = [ n.label for n in node.child_nodes() ]
            #logging.info(f'handling node {i}, parent of node {child1} and {child2}.')
            # calculate for each state the distance for this internal (parent) node
            # i: index of parent node
            # j: index of parent state
            # p: index of child node 1 state
            # q: index of child node 2 state
            # obtain the distance vector for the two child nodes
            distChild1 = distances[int(child1)]
            distChild2 = distances[int(child2)]
            for j,stateParent in enumerate(allStates):
                minDist, minP, minQ = float('inf'), 0, 0
                for p,stateChild1 in enumerate(allStates):
                    for q,stateChild2 in enumerate(allStates):
                        # distance incurred by the child1 branch
                        dist = abs(stateParent-stateChild1) + distChild1[p]
                        # distance incurred by the child2 branch
                        dist = dist + abs(stateParent-stateChild2) + distChild2[q]
                        if dist < minDist:
                            minDist, minP, minQ = dist, p, q
                # update parent distance/traceback vector
                distances[i][j] = minDist
                traceback[i][j] = [minP, minQ]
            '''print(distChild1)
            print(distChild2)
            print(distances[i])
            print(traceback[i])'''

    # final traceback
    bestStateIndex = {}
    # find index of the best state of the root node (i.e. the global best)
    rootBestDist = min(distances[numNodeTotal-1])
    bestStateIndex[numNodeTotal-1] = distances[numNodeTotal-1].index(rootBestDist)

    # traceback by preorder (parent comes before child)
    for parent in tree.preorder_node_iter():
        if parent.is_leaf():
            pass
        else:
            parentBestStateIndex = bestStateIndex[parent.label]
            child1, child2 = [ n.label for n in parent.child_nodes() ]
            #print(f'best child of node {parent.label}', traceback[parent.label][parentBestStateIndex])
            bestStateIndex[child1], bestStateIndex[child2] = traceback[parent.label][parentBestStateIndex]

    # by level order traversal (BFS)
    # change node taxons to the best state (internal node) or the allele (leaf node) for output
    # change node labels for switch checking
    for i,node in enumerate(tree.level_order_node_iter()):
        if node.is_leaf():
            sampleIdx = str(node.taxon).replace("'","")
            node.taxon = dendropy.Taxon(label=f"{i}[{alleles[int(sampleIdx)]}]")
            node.label = f'{i}-{alleles[int(sampleIdx)]}'
        else:
            node.taxon = dendropy.Taxon(label=f"{i}[{allStates[bestStateIndex[node.label]]}]")
            node.label = f'{i}-{allStates[bestStateIndex[node.label]]}'

    # check each root-to-leaf path for recurrent state switches
    recurrents, checkedIDs = [], []
    for leaf in tree.leaf_node_iter():
        path, current_node = [], leaf
        while current_node:
            path.append(current_node.label)
            current_node = current_node.parent_node
        path = path[::-1] # Reverse the path to start from the root
        ids = [ l.split('-')[0] for l in path ]
        states = [ l.split('-')[1] for l in path ]
        for degree,id in enumerate(ids):
            if degree < 1: continue # only check node from degree > 2
            if id in checkedIDs: continue # skip already checked nodes
            # if this is a recurrent state
            if states[degree] in states[:degree-1] and states[degree] != states[degree-1]:
                recurrents.append(id)
            checkedIDs.append(id)
    nRecurrent = str(len(recurrents))

    # change labels for output
    for i,node in enumerate(tree.preorder_node_iter()):
        node.label = None

    #for node in tree.postorder_node_iter(): print(node)
    out.write('\t'.join([chr,start,end,cons,nRecurrent,tree.as_string("newick")]))

'''for n,stateIndex in bestState.items():
    print(n,allStates[stateIndex])'''

out.close()

logging.info('End of Program\n')

