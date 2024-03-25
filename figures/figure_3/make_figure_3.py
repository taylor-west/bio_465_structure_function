import os
import math
import re
import ast
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt



TARGET_UNIPROT_ID = 'P50933'
PATH_TO_CLUSTERS_DATA = os.path.join(os.getcwd(), "datafiles")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
CLUSTER_CSV_FILEPATH = '/content/clusters_K00016.csv'

def generate_figure_3(target_uniprot_id, cluster_csv_filepath):
    clusters_df = pd.read_csv(cluster_csv_filepath)
    clusters_df.dropna(inplace=True)
    subset_clusters_df = clusters_df[clusters_df['uniprot_id'] == target_uniprot_id]


if __name__ == "__main__":
    generate_figure_3(TARGET_UNIPROT_ID, CLUSTER_CSV_FILEPATH)