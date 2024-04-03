import pandas as pd
import requests

#######################################################
# TODO: combine this with the kegg_ortholog.py file
#######################################################

def find_uniprot_gene_collections(pathway_orthologs_filepath, temporary_results_filepath, num_valid_prots_threshold = 10, annotation_score_threshold=3.0):
    orthologs_df = pd.read_csv(pathway_orthologs_filepath)
    ortholog_ids = set(orthologs_df['ko_id'])
    filtered_orthologs_df = pd.DataFrame(columns=['ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
    for ortholog_id in ortholog_ids:
        ortholog_uniprots_df = orthologs_df[orthologs_df['ko_id'] == ortholog_id]
        if len(ortholog_uniprots_df) < num_valid_prots_threshold:
            continue
        uniprot_ids = list(ortholog_uniprots_df['uniprot_id'])
        url = 'https://rest.uniprot.org/uniprotkb/search?query='
        for i in range(len(uniprot_ids)):
            if i == len(uniprot_ids) - 1:
                # don't add another 'OR' to the query
                url += f'accession:{uniprot_ids[i]}'
            else:
                # add another 'OR' to the query
                url += f'accession:{uniprot_ids[i]}+OR+'
        url += '&fields=accession,gene_names,protein_name,organism_name,annotation_score,ft_binding&format=tsv'

        response = requests.get(url)

        if response.status_code == 200:
            with open (temporary_results_filepath, 'w') as file:
                file.write(response.text)
                file.close()
        else:
            print(f'failed to retrieve uniprot data for ortholog id {ortholog_id}')
            break

        entry_data_df = pd.read_csv(temporary_results_filepath, sep='\t')
        entry_data_df.dropna(inplace=True, subset=['Binding site'])
        entry_data_df = entry_data_df[entry_data_df['Annotation'] >= annotation_score_threshold]
        print(len(entry_data_df))
        if len(entry_data_df) >= num_valid_prots_threshold:
            filtered_ids = list(entry_data_df['Entry'])
            for filtered_id in filtered_ids:
                row = orthologs_df[orthologs_df['uniprot_id'] == filtered_id]
                filtered_orthologs_df = pd.concat([filtered_orthologs_df, row], ignore_index=True)

    return filtered_orthologs_df