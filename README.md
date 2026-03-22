# Y-HaploMatch

Y-HaploMatch is a small web app for comparing a user's Y-chromosome SNP file to ancient individuals from the Allen Ancient DNA Resource (AADR).

The user enters a Y haplogroup and uploads a genotype file. The app then:
- finds SNPs relevant to that haplogroup
- checks which of those SNPs appear in the uploaded file
- compares them to ancient individuals in the database
- shows matches in a table and on a map

This project is meant to be run locally on `localhost`.

## What You Need

If you are starting on a new computer, install these first:

1. `Git`
2. `Conda` (Miniconda, Anaconda, or Miniforge)
3. `Node.js` and `npm`

## Download The Project

Clone the repo:

```bash
git clone https://github.com/nadinepi/y-haplomatch.git
cd y-haplomatch
```

## Data Files

To run the app, the you need to install the following SQLite database file:

- `data/ydna.db`

You can download this from:
https://github.com/nadinepi/y-haplomatch/releases/tag/database

Unzip it and put it in here:

```text
data/ydna.db
```

## Set Up The Backend

Create the conda environment:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate popgen
```

## Set Up The Frontend

Install the frontend packages:

```bash
npm --prefix frontend install
```

## Run The App

Start the backend in one terminal:

```bash
./backend.sh
```

Start the frontend in another terminal:

```bash
./frontend.sh
```

Then open the frontend URL printed in the terminal. It is usually:

```text
http://localhost:5173
```

The backend usually starts on port `5174`, but the scripts can choose another free port if needed.

## How To Use It

1. Enter a Y haplogroup, for the demo enter `I1`
2. Upload a genotype file (`txt`, `csv`, `tsv`, or `raw`). The demo user file, Demo_I1_DNA.txt, is provided at demo_input/Demo_I1_DNA.txt
3. Click `Find Matches`
4. Explore the results in the table and on the map

There is also a 'Run demo user' button in the app, which will do the above for you in 1 click.

The output should look like:

You can also find this image as demo_output/demo_output.png in the repo.

## Rebuilding The Database (Optional)

Most users do not need to do this.

If you want to rebuild `data/ydna.db`, first download the files from:

[Project files folder](https://lunduniversityo365-my.sharepoint.com/personal/er2374el_lu_se/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fer2374el_lu_se%2FDocuments%2FStudents%2FBINP29%20-%20DNA%20Sequencing%20Informatics%20II%2FStudent%20projects%2FProjetcFiles&ga=1)

Download these folders:

- `AADR_54.1`
- `AncientYDNA`
- `TestUsers`

Then put them inside the `data/` folder.

After that, download and run PLINK on the AADR files so you create the Y chromosome `.raw` file.

Use this PLINK command:

```bash
plink \
  --bfile /home/inf-35-2025/BINP29/popgen_project/ProjetcFiles/AADR_54.1/v54.1_1240K_public \
  --chr 24 \
  --filter-males \
  --out /home/inf-35-2025/BINP29/popgen_project/data/plink_out/aadr_chrY \
  --recode A
```

This should create a `plink_out` folder inside `data/`, with this file:

```text
data/plink_out/aadr_chrY.raw
```

Make sure that file is there. Then run:

```bash
conda activate popgen
python scripts/build_profiles.py
```

This uses the tree files, SNP files, metadata files, and `aadr_chrY.raw` to create a new `data/ydna.db`.

## Notes

- The app is designed for local use, not for deployment to a public server, but it is the future goal.
- Matching is based on haplogroup-relevant Y SNPs and position-based comparison to ancient individuals.
- Some user inputs may be resolved using the 2016 reference and some using the 2019 reference, depending on what the app recognizes.
