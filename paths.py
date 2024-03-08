import os

CWD = os.getcwd()
PATH_TO_DATAFILES = os.path.join(CWD, "../datafiles")
PATH_TO_DATAFILES_FROM_MAIN = "datafiles" #check this one idk. This all might need to be refined...
PATH_TO_FIGURE1 = os.path.join(PATH_TO_DATAFILES, "figure1")
PATH_TO_MUSCLE_DATA = os.path.join(PATH_TO_DATAFILES, "muscle_data")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
PATH_TO_EVAL_FILES = os.path.join(PATH_TO_DATAFILES, "eval_files")
PATH_TO_ORTHOLOG_UNIPROTS = os.path.join(PATH_TO_DATAFILES, "ortholog_uniprots")
PATH_TO_KO_IDS = os.path.join(PATH_TO_ORTHOLOG_UNIPROTS, "ko_ids")
PATH_TO_PATHWAYS = os.path.join(PATH_TO_ORTHOLOG_UNIPROTS, "pathways")
PATH_TO_PDB_FILES = os.path.join(PATH_TO_DATAFILES, "pdb_files")
