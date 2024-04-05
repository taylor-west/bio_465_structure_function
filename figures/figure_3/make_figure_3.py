
import os
import sys
import math
import ast
import re
import requests
from itertools import chain
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


FIGURE_3_RESULTS_DIRECTORY_PATH = os.path.join(os.getcwd(), "figures", "figure_3", "figure_3_results")
COMMUNITY_COLOR_OPTIONS = ['orange','magenta','lightgreen','pink','lightblue','black','red','blue', "green", "brown"]


def generate_figure_3(target_uniprot_id, cluster_csv_filepath, add_notes_to_fig: False):
    clusters_df = pd.read_csv(cluster_csv_filepath)
    clusters_df.dropna(inplace=True)
    subset_clusters_df = clusters_df[clusters_df['uniprot_id'] == target_uniprot_id]

    # make the adj_dict
    adj_dict = make_adj_dict(subset_clusters_df)
    print(f"adj_dict: {adj_dict}")

    # make the communities
    G = nx.Graph(adj_dict)
    communities = make_communities(G)
    print(f"communities: {communities}")

    # assign colors to the communities
    current_color_index = 0
    community_colors = {}

    for community in communities:
        community_colors[community] = COMMUNITY_COLOR_OPTIONS[current_color_index]
        current_color_index = current_color_index+1

    # get known binding sites from UniProt
    known_binding_site_groups_dict = get_known_binding_site_groups(target_uniprot_id)
    known_binding_sites_flat = list(chain.from_iterable(list(known_binding_site_groups_dict.values()))) # flatten all values from expected_binding_site_groups_dict

    print(f"known_binding_site_groups: {known_binding_site_groups_dict}")
    print(f"known_binding_sites: {known_binding_sites_flat}")

    # add labels to known binding sites that were predicted
    labels = {}
    for node in G.nodes():
        if node in known_binding_sites_flat:
            #set the node name as the key and the label as its value
            labels[node] = node

    # adds notes to the graph
    if add_notes_to_fig.lower() == 'true':
      target_ko_id = get_ko_id_from_cluster_filepath(cluster_csv_filepath)
      organisms_list_filename = clusters_df.organisms_list.unique()[0]
      num_known_binding_sites_found = len(labels)
      num_known_binding_sites = len(known_binding_sites_flat)
      num_predicted_binding_sites = len(G.nodes)
      num_known_binding_sites_grouped_correctly = get_num_nodes_grouped_correctly(communities, known_binding_site_groups_dict)
      num_communities_found = len(communities)

      plt.figure()
      annotate_graph(plt.gca(), target_uniprot_id, target_ko_id, organisms_list_filename, num_communities_found, known_binding_site_groups_dict, num_known_binding_sites_found, num_known_binding_sites, num_predicted_binding_sites, num_known_binding_sites_grouped_correctly)


    # Choose a layout algorithm for positioning nodes
    pos = nx.spring_layout(G, k=0.15, iterations=20)  # You can choose any layout algorithm you prefer

    # Draw the graph with colors based on community
    nx.draw(
        G,
        pos,
        with_labels=False,
        node_size=150,
        node_color=[get_color_for_node(node, communities, community_colors) for node in G.nodes()]
    )

    # Add labels to nodes in binding_site
    nx.draw_networkx_labels(G, pos, labels, font_size=10)

    # save the figure to a png
    save_file_name = "figure3_" + target_uniprot_id + ".png"
    save_file_path = os.path.join(FIGURE_3_RESULTS_DIRECTORY_PATH, save_file_name)
    plt.savefig(save_file_path)

    return


def make_adj_dict(subset_clusters_df):
  adj_dict = {}
  for i, row in subset_clusters_df.iterrows():
      cluster = ast.literal_eval(row['clusters'])
      adj_dict[row['seq_pos']] = cluster
  return adj_dict

def make_communities(G):
  return list(nx.community.greedy_modularity_communities(G))

# get known binding sites from UniProt
def get_known_binding_site_groups(uniprot_id):
  url = f'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession={uniprot_id}'

  response = requests.get(url, headers={ "Accept" : "application/json"})

  if response.status_code == 200:
      with open ('results.json', 'w') as file:
          file.write(response.text)
          file.close()
  else:
      print(f'failed to retrieve uniprot data for ortholog id {uniprot_id}')


  entry_data_df = pd.read_json('results.json')

  site_dict = {}

  for i in entry_data_df['features'][0]:
    if i['category'] == 'DOMAINS_AND_SITES':
      if i['type'] == 'BINDING':
        if i['ligand']['name'] not in site_dict.keys():
          site_dict[i['ligand']['name']] = [i]
        else:
          site_dict[i['ligand']['name']] += [i]

      # uncomment if you also want to look for active sites
      # else:
      #   if i['type'] not in site_dict.keys():
      #     site_dict[i['type']] = [i]
      #   else:
      #     site_dict[i['type']] += [i]

  known_binding_sites_dict = {}
  
  for ligand in site_dict:
    for residue in site_dict[ligand]:
      if ligand not in known_binding_sites_dict.keys():
        known_binding_sites_dict[ligand] = list(range(int(residue['begin']),int(residue['end'])+1))
      else:
        known_binding_sites_dict[ligand] += list(range(int(residue['begin']),int(residue['end'])+1))

  return known_binding_sites_dict

# assignes a color to each node based on the node's community
def get_color_for_node(target_node, communities, community_colors_dict):
  for community in communities:
    if target_node in community:
      return community_colors_dict[community]

# parses the KO_ID from the cluster filepath
def get_ko_id_from_cluster_filepath(cluster_csv_filepath):
  return re.findall('.*_(K.*).csv', cluster_csv_filepath)[0]

# returns the sum of the number of nodes correctly grouped into the community with the largest number of nodes belonging to each binding site
def get_num_nodes_grouped_correctly(communities, known_binding_site_groups):
  # loop through known binding site groups and count the largest group of that site in a single community
  binding_sites_biggest_correct_grouping = {}
  
  for binding_site_group_name in known_binding_site_groups.keys():
    curr_binding_site_group = known_binding_site_groups[binding_site_group_name]
    print(f"binding_site_group: {binding_site_group_name}  ({len(curr_binding_site_group)} nodes total)") # TODO: remove

    community_node_grouping_count_dict = {}
    for community in communities:
      community_node_grouping_count_dict[community] = 0

    print(f"  community_node_grouping_count_dict: {community_node_grouping_count_dict}") # TODO: remove

    # loop through nodes in the binding_site_group
    for node in curr_binding_site_group:

      # loop through communities and find nodes that are both known and predicted (same as labels)
      for community in communities:
        if node in community:
          community_node_grouping_count_dict[community] += 1 #community_node_grouping_count_dict[community] + 1
    
    # finds the community with the largest number of correctly grouped nodes for this binding site groups
    largest_correctly_grouped_community_name = max(community_node_grouping_count_dict, key=community_node_grouping_count_dict.get)
    num_largest_correctly_grouped = community_node_grouping_count_dict[largest_correctly_grouped_community_name]

    
    print(f"  winning community: {largest_correctly_grouped_community_name}") # TODO: remove
    print(f"  winning count: {num_largest_correctly_grouped}") # TODO: remove

    # sets the winning count for this binding site
    binding_sites_biggest_correct_grouping[binding_site_group_name] = num_largest_correctly_grouped
  
  total_correct_sum = sum(binding_sites_biggest_correct_grouping.values())
  print(f"binding_sites_biggest_correct_grouping = {binding_sites_biggest_correct_grouping}") # TODO: remove
  print(f"total_correct_sum = {total_correct_sum}") # TODO: remove
  return total_correct_sum
   

def annotate_graph(figure, target_uniprot_id, target_ko_id, organisms_list_filename, num_communities_found, known_binding_site_groups_dict, num_known_binding_sites_found, num_known_binding_sites, num_predicted_binding_sites, num_known_binding_sites_grouped_correctly):
  # top left corner
  figure.text(0.0, 1.0, f"uniprot_id={target_uniprot_id}",
      transform=plt.gca().transAxes)

  # top center
  figure.text(0.75, 1.0, f"ko_id={target_ko_id}",
      transform=plt.gca().transAxes)

  # bottom left corner
  figure.text(0.0, 0.00, f"organisms_list='{organisms_list_filename}'",
    transform=plt.gca().transAxes)
  
  figure.text(0.0, -0.05, f"known_binding_sites{known_binding_site_groups_dict}",
      transform=plt.gca().transAxes)

  figure.text(0.0, -0.125, f"{num_communities_found} communities found ({len(known_binding_site_groups_dict.keys())} binding site groups known)",
      transform=plt.gca().transAxes)
  
  percent_known_found = (round((num_known_binding_sites_found/num_known_binding_sites)*100,1) if num_known_binding_sites != 0 else "--")
  figure.text(0.0, -0.175, f"{num_known_binding_sites_found}/{num_known_binding_sites} ({percent_known_found}%) known sites found",
      transform=plt.gca().transAxes)
  
  percent_predicted_known = (round((num_known_binding_sites_found/num_predicted_binding_sites)*100,1) if num_predicted_binding_sites != 0 else "--")
  figure.text(0.0, -0.225, f"{num_known_binding_sites_found}/{num_predicted_binding_sites} ({percent_predicted_known}%) predicted sites are known",
      transform=plt.gca().transAxes)
  
  percent_known_grouped_correctly = (round((num_known_binding_sites_grouped_correctly/num_known_binding_sites_found)*100,1) if num_known_binding_sites_found != 0 else "--")
  figure.text(0.0, -0.275, f"{num_known_binding_sites_grouped_correctly}/{num_known_binding_sites_found} ({percent_known_grouped_correctly}%) known sites grouped correctly",
      transform=plt.gca().transAxes)
  return


if __name__ == "__main__":
    target_uniprot_id = sys.argv[1]
    cluster_csv_filepath = sys.argv[2]
    add_notes_to_fig = (sys.argv[3] if len(sys.argv) > 2 else None)

    generate_figure_3(target_uniprot_id, cluster_csv_filepath, add_notes_to_fig)
