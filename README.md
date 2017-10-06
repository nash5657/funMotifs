regMotifs is developed to annotate transcription factors (TF) with functional annotations. The main module is src/regMotifsMain.py
Parameters to specify the input files and other parameters should be given in a configuration file. An example of such configuration file is provided in conf/main_parameters.conf

The datafiles that are specified in the main_parameters.conf have to be available for the tool to run. Please follow ReadMe files listed in the ReadMe directory to generate input datafiles and annotations for running regMotifs.
A complete set of datafiles are prodvided to re-generate the annotated motifs that are reported in the current version of regMotifsDB. They can be downloaded on: bioinf.icm.uu.se/regMotifs or on ULAM (/data1/husen/datafiles.tar.gz)

The pipeline also creates a database and inserts the annotated motifs into a single table where each row represents a motif and the columns represent information about the motif (position, name, p-value, score) and the remaining columns are showing annotations per cell. In order for this task PostgreSQL has to be accessible on a host that has to be specified in the main parameters file (default is localhost). Also, database name and username parameters have to be specified. 

Python 2.7 with the following packages are required to run the pipeline:
- bedtools v25.9
- pybedtools (v0.7.8): for processing the annotation data files (e.g. conda install pybedtools=0.7.8 -c bioconda
- json: for processing dictionary files
- psycopg2 (v2.6.2): for connecting with the PostgreSQL database

The following packages are also required to run the helper modules:
- requests
- userlib 
- numpy 