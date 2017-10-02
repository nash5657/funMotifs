############################ Curate ENCODE Datasets (ChIP-seq and DNase-seq) #########################
Meatadata files are generated from the ENCODE project web portal. Each line in the file contains information about a data file including its format and accession code. The accession codes are used to download the files from https://www.encodeproject.org/files/.

#### Get ChIP-seq metadata file #####
ENCODE TF ChIP-seq datasets (copy the following URL and change the search features. After your customizations are done, copy the URL to get the path for a customized metadata file). e.g:
https://www.encodeproject.org/search/?searchTerm=chip-seq&status=released&assembly=hg19&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&target.investigated_as=transcription+factor&files.file_type=bed+narrowPeak&limit=all&frame=object&type=Experiment

In order to get the metadata file, append metadata.tsv to the URL and press enter from a browser or type wget "URL" in a terminal. e.g.:
wget "https://www.encodeproject.org/metadata/searchTerm=chip-seq&status=released&assembly=hg19&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&target.investigated_as=transcription+factor&files.file_type=bed+narrowPeak&limit=all&frame=object&type=Experiment/metadata.tsv"

mv metadata.tsv conf/metadataENCODETFChIPSeq29Nov2016.tsv
*Note '/' character should be replaced with '_' in the cell name column (#col 7); e.g :%s/NT2\/D1/NT2_D1/gc

The metadata file is parsed by ParseCellInfo.py to get the data files using the File Accession. The parser downloads the file through:
https://www.encodeproject.org/files/FileAccession/@@download/FileAccession

The metadata file path should be specified in conf/ParseCellInfo_params.conf file to inform the ParseCellInfo module. Assign the path to the ChIP-seq lines in the conf file.

#### Repeat the same steps to get DNase1 metadata file #####
#e.g:
https://www.encodeproject.org/batch_download/type%3DExperiment%26assay_slims%3DDNA%2Baccessibility%26assay_title%3DDNase-seq%26status%3Dreleased%26replicates.library.biosample.donor.organism.scientific_name%3DHomo%2Bsapiens%26assembly%3Dhg19%26files.file_type%3Dbed%2BnarrowPeak%26limit%3Dall
#Download the metadata file. e.g.:
wget "https://www.encodeproject.org/metadata/type=Experiment&assay_slims=DNA+accessibility&assay_title=DNase-seq&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=hg19&files.file_type=bed+narrowPeak&limit=all/metadata.tsv"

mv metadata.tsv metadataENCODEDNaseSeq29Nov2016.tsv

#Note '/' character should be replaced with '_' in the cell name column (#col 7); e.g :%s/NT2\/D1/NT2_D1/gc

The metadata file path should be specified in conf/ParseCellInfo_params.conf file to inform the ParseCellInfo module. Assign the path to the DNase-seq lines in the conf file.
