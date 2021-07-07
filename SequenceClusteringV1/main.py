import numpy as np
from KMedoid import *
from sklearn import metrics
import matplotlib.cm as cm
import matplotlib.pyplot as plt

def ParseFileToSortedTuples(filename):
    with open(filename, 'r') as reader:
        s = set(reader.readlines()) # remove duplicates
        return list(ParseStrings(s))
    
# Get a set of sequences in list form
def ParseStrings(string_set):
    list_set = set()
    for i in string_set:
        i = i.replace('\n', '')
        i = i.split(",")
        list_set.add(tuple(i)) # remove duplicates sequence
    return list_set


def silhouetteAnalysis(inputData):
    n_clusters_range = range(10, 41)
    for n_cluster in n_clusters_range:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)
        plt.suptitle(("Silhouette analysis for Kmedoid clustering on sample data "
                  "with n_clusters = %d" % n_cluster),
                 fontsize=14, fontweight='bold')

        # graph 1 is silhouette plot
        # it sets the x-axis range which corresponds to the silhouette coefficient
        ax1.set_xlim([-1, 1])
        # it sets the y-axis range which is the space for all clusters points 
        # 10*(n_cluster + 1) adds blank space between clusters
        ax1.set_ylim([0, len(inputData) + 10*(n_cluster + 1)])
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The sample's silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
        

        kme = Kmedoid(inputData, n_cluster)
        kme.cluster()
        clusters_labels = kme.labels
        dist_matrix     = kme.dist_matrix

        cluster_score = getQualityScore(dist_matrix, clusters_labels)
        print(cluster_score)
        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=cluster_score, color="red", linestyle="--")
        ax1.set_yticks([])  # Clear the yaxis labels / ticks

        # Compute the silhouette scores for each sample
        y_lower = 10
        sample_silhouette_values = metrics.silhouette_samples(dist_matrix, clusters_labels)
        for i in range(n_cluster):
            ith_cluster_sample_silhouette_values = sample_silhouette_values[clusters_labels == i]
            ith_cluster_sample_silhouette_values.sort()
            size_cluster_i = ith_cluster_sample_silhouette_values.shape[0]
            y_upper = size_cluster_i + y_lower
            color = cm.nipy_spectral(float(i) / n_cluster) # create a unique color for each cluster
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_sample_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples
    plt.show()





def getQualityScore(dist_matrix, labels):
    quality_score = metrics.silhouette_score(dist_matrix, labels, metric="precomputed")
    return quality_score



if __name__ == '__main__':
    inputData = ParseFileToSortedTuples("toy_vector_list.txt")
    # 1) print set of sequences
    print("total sequences: ", len(inputData))

    # 2) find best cluster_number 
    silhouetteAnalysis(inputData)
  
    

    ####################################
    # TO DO!
    # X continue coding __compute_medoids method
    # X test kMedoid and get it working with a score
    # x hyperparameter find better # of clusters
    # **use graph idea** 

    # Notes:
    # critique on find approximate intermediary seq for calculate centroid-seq distance
    # how to find centroids? - maybe not possible for string disimilarity, typically use euclidean distance function and Means of data points