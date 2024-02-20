# bio_465_structure_function
BIO 465 - Capstone Project


## Setup
- Install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) if you don't already have it
- run `conda create --name <env> --file requirements.txt && conda activate <env>`


### If not using arm64 architecture
- create a new conda environement by running `conda create --name <env>`
- activate this environment with `conda activate <env>`
- install the following packages:
   - `conda install python`
   - `conda install requests`
   - `conda install numpy`
   - `conda install pandas`
   - `conda install -c etetoolkit ete3`
   - `conda install biopython`

[post about conditional requirements files](https://stackoverflow.com/questions/29222269/is-there-a-way-to-have-a-conditional-requirements-txt-file-for-my-python-applica)



## Guidelines
### Branches
The `main` branch is protected to prevent direct pushes and commits. All work should be done on a separate branch and merged in via Pull Request. For style purposes, it is reccomended that you preface your branch with your name.

#### Workflow
1. Pull the latest changes from main
   - `git checkout main && git pull`
2. Create a new branch based off of main.
   - `git checkout -b <your_branch_name>`
   - For style purposes, it is reccomended that you preface your branch with your name (e.g.`user1/test_branch`).
3. Make and commit your changes.
4. Re-merge master and handle merge conflicts
   - `git checkout main && git pull && git checkout <your_branch_name> && git merge main`
5. Commit and push your changes to your branch
   - `git push --set-upstream origin <your_branch_name>`
6. [Create a Pull Request](https://github.com/taylor-west/bio_465_structure_function/pulls) on GitHub and write a brief description of your changes
7. Merge branch into `main`


