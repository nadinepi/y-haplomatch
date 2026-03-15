# Y HaploMatch

A web app for exploring Y-chromosome SNP matches with ancient DNA samples.

## Features
- upload your Y-chromosome genotype file
- see matches with ancient individuals from the AADR database
- interactive map and results table
- filter and explore shared SNPs

## Software Stack
- backend: Flask (Python), SQLite
- frontend: React, Vite, Leaflet
- see `backend.sh` and `frontend.sh` for running locally

## Usage
1. enter your Y haplogroup and upload your genotype file
2. click "Find Matches"
3. view results on the map and in the table


## Installation 

Please download 3 folders:
- AADR_54.1
- AncientYDNA
- TestUsers

from:
https://lunduniversityo365-my.sharepoint.com/:f:/r/personal/er2374el_lu_se/Documents/Students/BINP29%20-%20DNA%20Sequencing%20Informatics%20II/Student%20projects/data?csf=1&web=1&e=MsrxoD

And place them into the data/ folder. Please open the AADR Annotations 2025.xlsx and export it to a csv format.
Unzip files, and delete unneccessary ones until your structure looks like:

├── AADR_54.1
│   ├── AADR Annotations 2025.csv
│   ├── Ancient_samples.txt
│   ├── Modern_samples.txt
│   ├── v54.1_1240K_public.bed
│   ├── v54.1_1240K_public.bim
│   └── v54.1_1240K_public.fam
├── AncientYDNA
│   ├── chrY_hGrpTree_isogg2016.txt
│   ├── chrY_locusFile_b37_isogg2016.txt
│   └── snpFile_b37_isogg2019.txt
├── TestUsers
│   ├── Test1_DNA.txt
│   ├── Test1.txt
│   ├── Test2_DNA.csv
│   ├── Test2.txt
│   ├── Test3.csv
│   ├── Test3.txt
│   ├── Test4_DNA.txt
│   ├── Test4.txt
│   ├── Test5_DNA.txt
│   ├── Test5.txt

---

> this project is under active development. check back for updates!