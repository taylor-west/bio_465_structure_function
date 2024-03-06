import pandas as pd
import os
import cluster_analysis_3d.CalculateResidueDistanceWithDataframeInput as cluster

PATH_TO_PARENT = os.path.dirname(os.getcwd())
PATH_TO_DATAFILES = os.path.join(PATH_TO_PARENT, "datafiles")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
PATH_TO_EVAL_FILES = os.path.join(PATH_TO_DATAFILES, "eval_files")


def read_uniprot_file_to_analyze_active_sites(directory, filename):
    important_positions = {}
    act_str = 'ACT_SITE'
    bind_str = 'BINDING'
    important_positions[act_str] = []
    important_positions[bind_str] = []
    for file in os.listdir(directory):
        if file == filename:
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as inF:
                for line in file:
                    line = line.strip()
                    if line.startswith('FT'):
                        entry = line.split()
                        if entry[1] == act_str:
                            important_positions[act_str].append(entry[2])
                        if entry[1] == bind_str:
                            important_positions[bind_str].append(entry[2])
    return important_positions

def get_key_rows_only(dataframe):
    grouped_df = dataframe.groupby('entry_id')
    aggregated_df = pd.DataFrame()
    for entry_id, group in grouped_df:
        entry_id = group['entry_id'].iloc[0]
        filename = entry_id + ".txt"
        key_positions = read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, filename)
        only_key_rows_df = group[group['seq_pos'].isin(key_positions)]
        # Store key_positions in a list
        key_positions_list = [key_positions] * len(only_key_rows_df)
        # Add a new column 'key_positions' to the filtered DataFrame
        only_key_rows_df['key_positions'] = key_positions_list
        aggregated_df = pd.concat([aggregated_df, only_key_rows_df], ignore_index=True)
    return aggregated_df

def analyze_key_locations(dataframe):
    key_results_df = get_key_rows_only(dataframe)
    key_results_df['common_positions'] = None
    key_results_df['num_key_present'] = 0
    key_results_df['percent_key_present'] = None
    for index, row in key_results_df.iterrows():
        cluster_positions = set(row['clusters'])
        key_positions = set(row['key_positions'])
        common_positions = list(cluster_positions.intersection(key_positions))
        key_results_df.at[index, 'common_positions'] = common_positions
        num_key_present = len(common_positions)
        percent_key_present = len(common_positions) / len(cluster_positions)
        key_results_df.at[index, 'num_key_present'] = num_key_present
        key_results_df.at[index, 'percent_key_present'] = percent_key_present
    return key_results_df

def evaluate_results(dataframe):
    key_df = get_key_rows_only(dataframe)
    results_df = analyze_key_locations(key_df)
    results_df.to_csv(os.path.join(PATH_TO_EVAL_FILES, 'analyzed_clusters.csv'))

invariant_res_df = pd.read_csv('../datafiles/MSA_results.csv')
clusters_df = pd.read_csv('../datafiles/clusters.csv')
invariant_res_df.drop(columns='Unnamed: 0', inplace=True)

result = cluster.find_clusters(invariant_res_df, 10)
print(result)
result2 = cluster.filter_interesting_clusters(clusters_df, 20)
print(result2)
evaluate_results(result)







