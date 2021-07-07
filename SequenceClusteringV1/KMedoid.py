import random
import numpy as np
import sys
import sys
import math

#k-medoids minimizes the sum of dissimilarities between points labeled to be in a cluster and a point designated as the center of that cluster. 
#In contrast to the k-means algorithm, k-medoids chooses datapoints as centers ( medoids or exemplars).
class Kmedoid:
    def __init__(self, input, n_clusters, max_iter=100, random_state=123):
        self.input = input
        self.dist_matrix = self.ComputeDistanceMatrix()
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.random_state = random_state
        self.labels = np.zeros(len(self.input))
        self.medoids = self.__initialize_medoids()

    # get initial medoids by randomly select sequences from the input set 
    def __initialize_medoids(self):
        medoids = list()
        sequencesCount = len(self.input)
        randomIndices = set()
        random.seed(self.random_state)
        
        for i in range(0, self.n_clusters):
            foundIndice = False
            while(foundIndice == False):
                index = random.randrange(sequencesCount)
                if index not in randomIndices:
                    foundIndice = True
                    randomIndices.add(index)
                    medoids.append(self.input[index])
        #print("initial medoids: ", medoids)
        return medoids

    def cluster(self):
        self.__assign_seq_to_medoid()

        # iteratively calculate the medoids and re-assign medoids till the medoid stops changing
        last_medoids = list()
        iteration = 0
        while self.__IsmedoidsChanged(last_medoids) and iteration <= self.max_iter:
            last_medoids = self.medoids.copy()
            self.__compute_medoids()
            self.__assign_seq_to_medoid()
            iteration += 1
            
    # check if the medoids has changed compared to the last iteration's medoids
    def __IsmedoidsChanged(self, lastmedoids):
        isChanged = False
        # Not changed if the last medoids can be found the current medoids, 
        # and current medoids can be found in last medoids
        for last_medoid in lastmedoids:
            if(last_medoid not in self.medoids):
                isChanged = True
                break 
        for medoid in self.medoids:
            if(medoid not in lastmedoids):
                isChanged = True
                break 
        return isChanged

    # populate labels attribute
    def __assign_seq_to_medoid(self):
        # create and populate a distance matrix of sequences vs medoids, of shape( number_sequences, number_medoids)
        medoid_distanceMatrix = np.zeros( (len(self.input), self.n_clusters) )
        for medoid in self.medoids:
            self.__computemedoidDistance(medoid, medoid_distanceMatrix)

        # assign labels for sequences
        seq_index = 0
        for row in medoid_distanceMatrix:
            min_dist = min(row)
            indexes = np.where(row == min_dist)
            self.labels[seq_index]=indexes[0][0]
            seq_index += 1
        
    def __computemedoidDistance(self, medoid, medoid_distanceMatrix):
        medoid_medoids_index = self.medoids.index(medoid)
        medoid_input_index   = self.input.index(medoid)

        for seq_i in range(0, len(self.input)):
            # the distance between a sequence and medoid
            medoid_distanceMatrix[seq_i, medoid_medoids_index] = self.dist_matrix[seq_i, medoid_input_index]
    
    # Generate a distance matrix between each sequences - An array of pairwise distances between sequences, array-like of shape (n_samples_a, n_samples_a)
    # the distance is a edit distance via Longest common substring method
    def ComputeDistanceMatrix(self):
        sequencesCount = len(self.input)
        distanceMatrix = np.zeros( (sequencesCount, sequencesCount) )
        # only compute the distance between sequence i and j once aka: (j,i) and (i,j), dstiance between sequence i and i is default to 0
        for i in range(0, sequencesCount):
            for j in range(0, sequencesCount):
                if(self.input[i] != self.input[j] and distanceMatrix[j][i] == 0):
                    dist = self.__calculate_LCS_EditDistance(self.input[i], self.input[j])
                    distanceMatrix[j][i] = dist
                    distanceMatrix[i][j] = dist
        return distanceMatrix

    def __calculate_LCS_EditDistance(self, curSeq, otherSeq):
            otherNodeLen = len(otherSeq)
            curNodeLen   = len(curSeq)
            commonCount = 0
            currentNodeQueue = list(curSeq)
            for ele in otherSeq:
                for i in range(len(currentNodeQueue)):
                    if(ele == currentNodeQueue[i]):
                        commonCount += 1
                        del currentNodeQueue[:i+1]
                        break
            return (otherNodeLen - commonCount) + (curNodeLen - commonCount)

    def __compute_medoids(self):
        # within each cluster, for each seq, sum the distances from all other elements to itself, aka:
            #   1    2    3
            # 1 0.0  2.0  1.0
            # 2 2.0  0.0  2.0
            # 3 1.0  2.0  0.0
            #------------------
           #sum 3.0  4.0  3.0

        # new medoids' indices
        medoids_indices = list()
        # create a sequence vs distance sum matrix for each cluster
        for i in range(0, self.n_clusters):
            cluster_seq_indices = np.where(self.labels == i) # find indices of sequences belonging to the ith cluster
            seq_count = len(cluster_seq_indices[0])
            cluster_seq_dist_sum = np.zeros( (seq_count, 2) ) # create a 1D array to store the sum for each sequence in the cluster
            # sum distance for each sequence 
            for j in range(0, seq_count):
                indice = cluster_seq_indices[0][j]
                sum = 0
                for other_i in cluster_seq_indices[0]:
                    sum += self.dist_matrix[other_i, indice]
                cluster_seq_dist_sum[j][0] = indice
                cluster_seq_dist_sum[j][1] = sum
            
            # find medoid whose distances is smallest
            #print("A cluster's distance sums for each sequence: ", cluster_seq_dist_sum)
            min_sum = min(cluster_seq_dist_sum[:,1])
            min_sum_index = np.where(cluster_seq_dist_sum[:,1] == min_sum)
            new_medoid_index = cluster_seq_dist_sum[min_sum_index[0][0],0].astype(int)
            medoids_indices.append(new_medoid_index)
        # set medoids based on indices
        self.medoids.clear()
        for medoid_index in medoids_indices:
            self.medoids.append(self.input[medoid_index])
        #print("new Medoids: ", self.medoids)




            















    

    #def cluster(self, graph, medoids):
    #    # calculate the distances between nodes in graph and medoids, stored in distance matrix
    #    nodes = graph.getNodes()
    #    distanceMatrix = np.zeros( (len(nodes), len(medoids)) )
    #    for i in range(0, len(nodes)):
    #        for j in range(0, len(medoids)):
    #            distance = graph.calculateEditDistance(nodes[i], medoids[j])
    #            distanceMatrix[i][j] = distance 
        
    #    # find the medoid for each node by smallest disimilarity
    #    node_index = 0
    #    for row in distanceMatrix:
    #        minDist = sys.maxsize
    #        minDistIndex = 0
    #        index = 0
    #        for dist in row:
    #            if dist < minDist:
    #                minDist = dist
    #                minDistIndex = index
    #            index +=1
    #        medoid = medoids[minDistIndex]
    #        nodes[node_index].addAdjacentNode(medoid, minDist)
    #        node_index += 1



