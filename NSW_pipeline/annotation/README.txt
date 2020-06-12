################################################################################
README for ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
           ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/
           ftp://ftp.ncbi.nlm.nih.gov/genomes/all/

Last updated: November 01, 2019
################################################################################

==========
Background
==========
Sequence data is provided for all single organism genome assemblies that are 
included in NCBI's Assembly resource (www.ncbi.nlm.nih.gov/assembly/).  This 
includes submissions to databases of the International Nucleotide Sequence 
Database Collaboration, which are available in NCBI's GenBank database, as well 
as the subset of those submissions that are included in NCBI's RefSeq Genomes 
project. 

Available by anonymous FTP at:
     ftp://ftp.ncbi.nlm.nih.gov/genomes/

Please refer to README files and the FTP FAQ for additional information:
     https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/

Subscribe to the genomes-announce mail list to be informed of changes to the
NCBI genomes FTP site:
     https://www.ncbi.nlm.nih.gov/mailman/listinfo/genomes-announce


=====================================================================
Genome sequence and annotation data is provided in three directories:
=====================================================================
1) all:     content is the union of GenBank and RefSeq assemblies. 
            Two directories under "all" are named for the accession prefix (GCA
            or GCF) and these directories contain another three levels of 
            directories named for digits 1-3, 4-6 & 7-9 of the assembly 
            accession. The next level is the data directories for individual 
            assembly versions. Only data directories for "latest" assemblies
            are refreshed when annotation is updated or when software updates
            are released, so new file formats or improvements to existing 
            formats are not available for non-latest assemblies.
            A third directory, named "annotation_releases" contains the products
            of the NCBI Eukaryotic Genome Annotation Pipeline (see below). The 
            data are organized first by taxonomy ID and then by annotation 
            release ID. It is expected that many users will prefer to access the
            annotation release data using the paths under the "refseq" directory
            that use the organism name (see below). 
2) genbank: content includes primary submissions of assembled genome sequence 
            and associated annotation data, if any, as exchanged among members 
            of the International Nucleotide Sequence Database Collaboration, 
            of which NCBI's GenBank database is a member. The GenBank directory 
            area includes genome sequence data for a larger number of organisms 
            than the RefSeq directory area; however, some assemblies are 
            unannotated. The sub-directory structure includes:
            a. archaea
            b. bacteria
            c. fungi
            d. invertebrate
            e. metagenomes
            f. other -  this directory includes synthetic genomes
            g. plant
            h. protozoa
            i. vertebrate_mammalian
            j. vertebrate_other
            k. viral
3) refseq:  content includes assembled genome sequence and RefSeq annotation 
            data. All prokaryotic and eukaryotic RefSeq genomes have annotation. 
            RefSeq annotation data may be calculated by NCBI annotation  
            pipelines or propagated from the GenBank submission. The RefSeq 
            directory area includes fewer organisms than the GenBank directory
            area because not all genome assemblies are selected for the RefSeq
            project.
            Sub-directories include:
            a. archaea
            b. bacteria
            c. fungi
            d. invertebrate
            e. plant
            f. protozoa
            g. vertebrate_mammalian
            h. vertebrate_other 
            i. viral
            j. mitochondrion [Content of the mitochondrion, plasmid and plastid
            k. plasmid     directories is from the RefSeq release FTP site. See 
            l. plastid     ftp://ftp.ncbi.nlm.nih.gov/refseq/release/README]

Data are further organized within each of the above directories as a series of 
directories named as the species binomial. For example:
   ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/
           - or - 
   ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/

The next hierarchy provides access to all assemblies for the species, latest 
assemblies, and selected reference or representative assemblies for the species 
(if any). Within these groupings, sequence and annotation (and other) data is 
provided per assembly in a series of directories that are named using the rule:

   [Assembly accession.version]_[assembly name]

For example, the directory hierarchy for the GenBank Bacillus thuringiensis 
strain 97-27 genome, which has the assembly accession GCA_000008505.1 and 
default assembly name ASM850v1 looks like this:  
   /genomes/genbank/bacteria/Bacillus_thuringiensis/all_assembly_versions/GCA_000008505.1_ASM850v1

The directory hierarchy for the RefSeq annotated human reference genome which 
has the assembly accession GCF_000001405.39 and assembly name GRCh38.p13 looks 
like this:
   /genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13

Species that have been annotated by the NCBI Eukaryotic Genome Annotation 
Pipeline will also have a directory named "annotation_releases" described below.

Genome assemblies of interest can be identified using the NCBI Assembly resource
(www.ncbi.nlm.nih.gov/assembly), or by using the assembly summary report files 
that are provided for both all genbank and all refseq assemblies:
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
or ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
or ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

Assembly summary report files containing information on assemblies for a 
particular taxonomic group or species are provided in the group and 
Genus_species directories under the "genbank" and "refseq" directory trees. e.g.
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/Sulfolobus_islandicus/assembly_summary.txt

Search the meta-data fields, or filter the files, to find assemblies of 
interest.


===========================
Data provided per assembly:
===========================
Sequence and other data files provided per assembly are named according to the 
rule:
[assembly accession.version]_[assembly name]_[content type].[optional format]

File formats and content:

   assembly_status.txt
       A text file reporting the current status of the version of the assembly
       for which data is provided. Any assembly anomalies are also reported.
   *_assembly_report.txt file
       Tab-delimited text file reporting the name, role and sequence 
       accession.version for objects in the assembly. The file header contains 
       meta-data for the assembly including: assembly name, assembly 
       accession.version, scientific name of the organism and its taxonomy ID, 
       assembly submitter, and sequence release date.
   *_assembly_stats.txt file
       Tab-delimited text file reporting statistics for the assembly including: 
       total length, ungapped length, contig & scaffold counts, contig-N50, 
       scaffold-L50, scaffold-N50, scaffold-N75, and scaffold-N90
   *_assembly_regions.txt
       Provided for assemblies that include alternate or patch assembly units. 
       Tab-delimited text file reporting the location of genomic regions and 
       listing the alt/patch scaffolds placed within those regions.
   *_assembly_structure directory
       This directory will only be present if the assembly has internal 
       structure. When present, it will contain AGP files that define how 
       component sequences are organized into scaffolds and/or chromosomes. 
       Other files define how scaffolds and chromosomes are organized into 
       non-nuclear and other assembly-units, and how any alternate or patch 
       scaffolds are placed relative to the chromosomes. Refer to the README.txt
       file in the assembly_structure directory for additional information.
   *_cds_from_genomic.fna.gz
       FASTA format of the nucleotide sequences corresponding to all CDS 
       features annotated on the assembly, based on the genome sequence. See 
       the "Description of files" section below for details of the file format.
   *_feature_count.txt.gz
       Tab-delimited text file reporting counts of gene, RNA, CDS, and similar
       features, based on data reported in the *_feature_table.txt.gz file.
       See the "Description of files" section below for details of the file 
       format.
   *_feature_table.txt.gz
       Tab-delimited text file reporting locations and attributes for a subset 
       of annotated features. Included feature types are: gene, CDS, RNA (all 
       types), operon, C/V/N/S_region, and V/D/J_segment. Replaces the .ptt & 
       .rnt format files that were provided in the old genomes FTP directories.
       See the "Description of files" section below for details of the file 
       format.
   *_genomic.fna.gz file
       FASTA format of the genomic sequence(s) in the assembly. Repetitive 
       sequences in eukaryotes are masked to lower-case (see below).
       The FASTA title is formatted as sequence accession.version plus 
       description. The genomic.fna.gz file includes all top-level sequences in
       the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
       unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
       that are part of the chromosomes are not included because they are
       redundant with the chromosome sequences; sequences for these placed 
       scaffolds are provided under the assembly_structure directory.
   *_genomic.gbff.gz file
       GenBank flat file format of the genomic sequence(s) in the assembly. This
       file includes both the genomic sequence and the CONTIG description (for 
       CON records), hence, it replaces both the .gbk & .gbs format files that 
       were provided in the old genomes FTP directories.
   *_genomic.gff.gz file
       Annotation of the genomic sequence(s) in Generic Feature Format Version 3
       (GFF3). Sequence identifiers are provided as accession.version.
       Additional information about NCBI's GFF files is available at 
       ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.
   *_genomic.gtf.gz file
       Annotation of the genomic sequence(s) in Gene Transfer Format Version 2.2
       (GTF2.2). Sequence identifiers are provided as accession.version.
   *_genomic_gaps.txt.gz
       Tab-delimited text file reporting the coordinates of all gaps in the 
       top-level genomic sequences. The gaps reported include gaps specified in
       the AGP files, gaps annotated on the component sequences, and any other 
       run of 10 or more Ns in the sequences. See the "Description of files" 
       section below for details of the file format.
   *_protein.faa.gz file
       FASTA format sequences of the accessioned protein products annotated on
       the genome assembly. The FASTA title is formatted as sequence 
       accession.version plus description.
   *_protein.gpff.gz file
       GenPept format of the accessioned protein products annotated on the 
       genome assembly
   *_rm.out.gz file
       RepeatMasker output; 
       Provided for Eukaryotes 
   *_rm.run file
       Documentation of the RepeatMasker version, parameters, and library; 
       Provided for Eukaryotes 
   *_rna.fna.gz file
       FASTA format of accessioned RNA products annotated on the genome 
       assembly; Provided for RefSeq assemblies as relevant (Note, RNA and mRNA 
       products are not instantiated as a separate accessioned record in GenBank
       but are provided for some RefSeq genomes, most notably the eukaryotes.)
       The FASTA title is provided as sequence accession.version plus 
       description.
   *_rna.gbff.gz file
       GenBank flat file format of RNA products annotated on the genome 
       assembly; Provided for RefSeq assemblies as relevant
   *_rna_from_genomic.fna.gz
       FASTA format of the nucleotide sequences corresponding to all RNA 
       features annotated on the assembly, based on the genome sequence. See 
       the "Description of files" section below for details of the file format.
   *_translated_cds.faa.gz
       FASTA sequences of individual CDS features annotated on the genomic 
       records, conceptually translated into protein sequence. The sequence 
       corresponds to the translation of the nucleotide sequence provided in the
       *_cds_from_genomic.fna.gz file. 
   *_wgsmaster.gbff.gz
       GenBank flat file format of the WGS master for the assembly (present only
       if a WGS master record exists for the sequences in the assembly).
   annotation_hashes.txt
       Tab-delimited text file reporting hash values for different aspects
       of the annotation data. See the "Description of files" section below 
       for details of the file format.
   md5checksums.txt file
       file checksums are provided for all data files in the directory

Additional directories and files provided for organisms annotated by the NCBI 
Eukaryotic Genome Annotation Pipeline:

   *_pseudo_without_product.fna.gz
       FASTA format of the genomic sequence corresponding to pseudogene and 
       other gene regions which do not have any associated transcribed RNA 
       products or translated protein products. It includes annotated gene 
       regions that require rearrangement to provide the final product, e.g.
       immunoglobulin segments. These sequences are not assigned accession 
       numbers, and are derived directly from the assembled genomic sequences.
       The FASTA title has a local sequence identifier, the Gene ID and gene 
       name.

Evidence_alignments directory
   *_cross_species_tx_alns.gff.gz
       Alignments of cDNAs, ESTs and TSAs from other species to the genomic
       sequence(s) in Generic Feature Format Version 3 (GFF3) [not all 
       annotation releases have cross-species alignments]. These alignments may
       have been used as evidence for gene prediction by the annotation 
       pipeline. Sequence identifiers are provided as accession.version. 
       Additional information about NCBI's GFF files is available at 
       ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.

   *_same_species_tx_alns.gff.gz
       Alignments of same-species cDNAs, ESTs and TSAs to the genomic 
       sequence(s) in Generic Feature Format Version 3 (GFF3). These alignments
       were used as evidence for gene prediction by the annotation pipeline. 
       Sequence identifiers are provided as accession.version. Additional 
       information about NCBI's GFF files is available at 
       ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.

Gnomon_models directory 
   *_gnomon_model.gff.gz
       Gnomon annotation of the genomic sequence(s) in Generic Feature Format
       Version 3 (GFF3). Sequence identifiers are provided as accession.version
       for the genomic sequences and Gnomon identifiers for the Gnomon models:
       gene.XXX for genes, GNOMON.XXX.m for transcripts and GNOMON.XXX.p for 
       proteins. These identifiers are NOT universally unique. They are unique
       per annotation release only. Additional information about NCBI's GFF 
       files is available at 
       ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.   
   *_gnomon_protein.faa.gz
       FASTA format sequences of Gnomon protein models annotated on the genome
       assembly. The FASTA title is the Gnomon identifier for the protein model
       (>gnl|GNOMON|XXX.p)   
   *_gnomon_rna.fna.gz
       FASTA format sequences of Gnomon transcript models annotated on the 
       genome assembly. The FASTA title is the Gnomon identifier for the 
       transcript (>gnl|GNOMON|XXX.m)

RefSeq_transcripts_alignments directory
   *_knownrefseq_alns.bam
       Alignments of the annotated Known RefSeq transcripts (identified with 
       accessions prefixed with NM_ and NR_) to the genome in BAM format [not 
       all annotation releases have Known RefSeq transcripts]. For more 
       information about the BAM format see: 
       https://samtools.github.io/hts-specs/SAMv1.pdf
   *_knownrefseq_alns.bam.bai
       Index of the BAM alignments of the annotated Known RefSeq transcripts 
       to the genome. [not all annotation releases have Known RefSeq 
       transcripts]
   *_modelrefseq_alns.bam
       Alignments of the annotated Model RefSeq transcripts (identified with 
       accessions prefixed with XM_ and XR_) to the genome in BAM format. For 
       more information about the BAM format see:
       https://samtools.github.io/hts-specs/SAMv1.pdf
   *_modelrefseq_alns.bam.bai
       Index of the BAM alignments of the annotated Model RefSeq transcripts to
       the genome.

Annotation_comparison directory
       This directory is only provided for re-annotations of the same species.
   *_compare_prev.txt.gz
       Matching genes and transcripts in the current and previous annotation 
       releases binned by type of difference (column 1 for genes and column 14 
       for transcripts), in tabular format.
   *_compare_prev.gbp.gz
       Genome Workbench project file for visualization and search of differences
       between the current and previous annotation releases.
       See how to download and use the 64-bit version of Genome Workbench:
       https://www.ncbi.nlm.nih.gov/tools/gbench/


=====================================
Data provided per annotation release:
=====================================
The annotation_releases directory offers data grouped by organism and specific 
annotation release (100, 101, etc.) for organisms that have been annotated 
by the NCBI Eukaryotic Genome Annotation Pipeline. 

Each annotation release corresponds to an annotation run. The annotation 
release identifiers (AR) are numbered sequentially starting at 100,
independently of the assembly used. An assembly may have been annotated multiple
times, and be featured in different annotation release directories. For example
Apis mellifera AR 103 was executed on the same assembly as A. mellifera AR 102,
Amel_4.5, using experimental evidence not available at the time AR 102 was 
produced. A. mellifera AR 104 was executed in 2018 on a newer assembly, 
Amel_HAv3.1.

The 'current' directory contains the data for the most recent annotation.
For many organisms, only the most recent annotation may be available. Previous 
annotations are available at 
ftp://ftp.ncbi.nlm.nih.gov/genomes/<organism>

For a small set of organisms including human (taxid 9606), we provide annotation
updates named <AR>.<date> that incorporate improvements made to genes and 
transcripts by RefSeq curation experts. See more details in:
https://ncbiinsights.ncbi.nlm.nih.gov/2019/03/26/human-genome-annotation-bimonthly-update/

Each annotation release directory contains:

README_[organism_name]_annotation_release_[annotation_release_id]
   This file provides information specific to the specific annotation release,
   including data freeze dates, release date and release number, and the 
   annotated assemblies.

[organism name]_ARXXX_annotation_report.xml
   This file is the XML version of the HTML report for the organism:
   https://www.ncbi.nlm.nih.gov/genome/annotation_euk/[org_name]/[annotation_release_id]/
   e.g. https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Homo_sapiens/108/
   It contains information on the annotation release, including:
      Important dates associated with the annotation
      Assemblies
      Gene and feature statistics
      Masking results
      Transcript and protein alignments used for the annotation
      Assembly-assembly alignments used to track genes from the previous 
      assembly to the current, or from the reference to an alternate assembly
      if relevant

One directory for each genome assembly that was annotated in the release. 
   Named as [assembly accession.version]_[assembly name]. 
   This directory contains the files provided for all genome assemblies plus 
   those files described above under "Additional directories and files provided
   for organisms annotated by the NCBI Eukaryotic Genome Annotation Pipeline".


=====================
Description of files:
=====================

Masking of fasta sequences in genomic.fna.gz files
--------------------------------------------------
Repetitive sequences in eukaryotic genome assembly sequence files, as 
identified by WindowMasker (Morgulis A, Gertz EM, Schaffer AA, Agarwala R. 
2006. Bioinformatics 22:134-41), have been masked to lower-case.

Alignment programs typically have parameters that control whether the program 
will ignore lower-case masking, treat it as soft-masking (i.e. only for finding 
initial matches) or treat it as hard-masking. By default NCBI BLAST will ignore 
lower-case masking but this can be changed by adding options to the blastn 
command-line.
To have blastn treat lower-case masking in the query sequence as soft-masking 
add:
     -lcase_masking
To have blastn treat lower-case masking in the query sequence as hard-masking 
add:
     -lcase_masking -soft_masking false

Alternatively, commands such as the following can be used to generate either 
unmasked sequence or sequence masked with Ns.

Example commands to remove lower-case masking:
perl -pe '/^[^>]/ and $_=uc' genomic.fna > genomic.unmasked.fna
  -or-
awk '{if(/^[^>]/)$0=toupper($0);print $0}' genomic.fna > genomic.unmasked.fna

Example commands to convert lower-case masking to masking with Ns (hard-masked):
perl -pe '/^[^>]/ and $_=~ s/[a-z]/N/g' genomic.fna > genomic.N-masked.fna
  -or-
awk '{if(/^[^>]/)gsub(/[a-z]/,"N");print $0}' genomic.fna > genomic.N-masked.fna


*_cds_from_genomic.fna.gz & *_rna_from_genomic.fna.gz
-----------------------------------------------------
FASTA sequences of individual features annotated on the genomic records. The 
sequences are based solely on the genome sequence and annotated feature at a
particular location. They may differ from the product sequences found in the 
*_rna.fna.gz and *_protein.faa.gz files which may be based on transcript or 
other data sources and include mismatches, indels, or additional sequence not 
found at a particular genomic location.

Seq-ids are constructed based on the following rule to ensure uniqueness:
lcl|<genomic accession.version>_<feature_type>_<product accession.version>_<counter>
Note the seq-id is not intended to be stable if the annotation is updated; in 
particular, addition or removal of feature(s) will cause the counter to change 
on following features.

The remainder of the FASTA definition line is composed of a series of qualifiers
bounded by brackets, as described at:
  https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html
  The qualifiers that may appear in these files are:
      gene
      locus_tag
      db_xref
      protein
      product
      ncRNA_class
      pseudo
      pseudogene
      frame
      partial
      transl_except
      exception
      protein_id
      location

Note that some qualifier values such as product names may themselves contain 
un-escaped brackets, which should be allowed for if parsing the files.
   
For CDS features that begin in frame 2 or 3, the first 1 or 2 bp of sequence
are trimmed from the CDS FASTA so that it always begins with the first complete
codon. The location and frame qualifiers are left unaltered; consequently, the 
length of the ranges in the location string may be 1-2 bp longer than the FASTA 
sequence.

For RefSeq assemblies annotated by NCBI's Eukaryotic Genome Annotation 
Pipeline, a gene may have a frameshifting indel(s) in the genome that is 
thought to result from a genome sequencing error; in these cases, the gene is 
still considered to be protein-coding and annotated with mRNA and CDS features, 
but the genome sequence won't translate correctly downstream from the 
frameshift. To compensate, the FASTA sequence of the genomic CDS and RNA 
features is modified with 1-2 bp gaps (aka "micro-introns") in order to 
restore the predicted reading frame. This modification is reflected by 1-2 bp 
micro-introns in the location qualifier. An equivalent modification is also
made in the *_genomic.gff.gz file. A protein-coding gene may also be annotated
with a CDS feature containing an in-frame stop codon that is translated as a
selenocysteine, subject to stop-codon readthrough, or thought to result from a
genome sequencing error; in these cases, a transl_except qualifier is provided
indicating the genomic location of the stop codon and its proposed translation.
For more details, see the section on "Annotation accommodations for putative 
assembly errors" in:
ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

Pseudogenes annotated with CDS features may be included in the 
*_cds_from_genomic.fna.gz file, and have FASTAs that are disrupted by 
frameshifting indels or in-frame stop codons. Pseudogene features can be
identified and screened out based on the presence of a [pseudo=true] qualifier
in the defline.


*_feature_count.txt.gz
----------------------
Tab-delimited text file reporting counts of gene, RNA, CDS, and similar 
features, based on data reported in the *_feature_table.txt.gz file (see below).
Separate counts are provided for different sets of sequences in the assembly 
corresponding to the primary assembly, non-nuclear assembly, all alt-loci 
sequences, and all patch scaffolds.

The file is tab delimited (including a #header) with the following columns:
col 1: Feature: INSDC feature type
col 2: Class: Gene features are subdivided into classes according to the gene 
       biotype. ncRNA features are subdivided according to the ncRNA_class. CDS 
       features are subdivided into with_protein and without_protein, depending 
       on whether the CDS feature has a protein accession assigned or not. CDS 
       features marked as without_protein include CDS features for C regions and 
       V/D/J segments of immunoglobulin and similar genes that undergo genomic 
       rearrangement, and pseudogenes.
col 3: Full Assembly: assembly accession.version for the full assembly
col 4: Assembly-unit accession: assembly accession.version for the assembly 
       unit.
col 5: Assembly-unit name: name of the assembly unit or set of sequences. For 
       assemblies with alt-loci or patch scaffolds, such as GRCh38.p11, all 
       sequences from all alt-loci or patches are combined together.
col 6: Unique Ids: counts of unique identifiers. For gene features, this is the
       count of unique GeneID db_xrefs, or locus_tags, such that genes that are
       annotated at more than one location on the assembly unit (e.g. on both 
       chrX and chrY in the PAR region) are counted once. For RNA and CDS 
       features, this is the count of unique product accessions. If no product 
       accession is assigned, such as for RNA features in GenBank genomes or CDS
       features classified as without_protein, then "na" is reported
col 7: Placements: count of all features of that type on the indicated assembly 
       unit or set of sequences.

Stats of common interest are:
- the count of protein-coding genes in the nuclear genome, which corresponds to
  "gene" in column 1, "protein_coding" in column 2, "Primary Assembly" in 
  column 5, and the count of Unique Ids as reported in column 6
- the count of distinct protein sequences annotated in the nuclear genome, which
  corresponds to "CDS" in column 1, "with_protein" in column 2, and the count of
  Unique Ids as reported in column 6
- the count of total CDS features with proteins annotated in the primary 
  assembly, regardless of whether two CDSes encode exactly the same protein and
  use the same RefSeq WP_ protein accession, which corresponds to "CDS" in 
  column 1, "with_protein" in column 2, and the count of Placements as reported
  in column 7


*_feature_table.txt.gz
----------------------
Tab-delimited text file reporting locations and attributes for a subset of 
annotated features. Included feature types are: gene, CDS, RNA (all types), 
operon, C/V/N/S_region, and V/D/J_segment. 

The file is tab delimited (including a #header) with the following columns:
col 1: feature: INSDC feature type
col 2: class: Gene features are subdivided into classes according to the gene 
       biotype computed based on the set of child features for that gene. See 
       the description of the gene_biotype attribute in the GFF3 documentation
       for more details: ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt
       ncRNA features are subdivided according to the ncRNA_class. CDS features
       are subdivided into with_protein and without_protein, depending on 
       whether the CDS feature has a protein accession assigned or not. CDS 
       features marked as without_protein include CDS features for C regions and 
       V/D/J segments of immunoglobulin and similar genes that undergo genomic 
       rearrangement, and pseudogenes.
col 3: assembly: assembly accession.version
col 4: assembly_unit: name of the assembly unit, such as "Primary Assembly", 
       "ALT_REF_LOCI_1", or "non-nuclear"
col 5: seq_type: sequence type, computed from the "Sequence-Role" and 
       "Assigned-Molecule-Location/Type" in the *_assembly_report.txt file. The
       value is computed as:
       if an assembled-molecule, then reports the location/type value. e.g. 
       chromosome, mitochondrion, or plasmid
       if an unlocalized-scaffold, then report "unlocalized scaffold on <type>".
       e.g. unlocalized scaffold on chromosome
       else the role, e.g. alternate scaffold, fix patch, or novel patch
col 6: chromosome
col 7: genomic_accession
col 8: start: feature start coordinate (base-1). start is always less than end
col 9: end: feature end coordinate (base-1)
col10: strand
col11: product_accession: accession.version of the product referenced by this 
       feature, if exists
col12: non-redundant_refseq: for bacteria and archaea assemblies, the 
       non-redundant WP_ protein accession corresponding to the CDS feature. May
       be the same as column 11, for RefSeq genomes annotated directly with WP_
       RefSeq proteins, or may be different, for genomes annotated with 
       genome-specific protein accessions (e.g. NP_ or YP_ RefSeq proteins) that
       reference a WP_ RefSeq accession.
col13: related_accession: for eukaryotic RefSeq annotations, the RefSeq protein
       accession corresponding to the transcript feature, or the RefSeq 
       transcript accession corresponding to the protein feature.
col14: name: For genes, this is the gene description or full name. For RNA, CDS,
       and some other features, this is the product name.
col15: symbol: gene symbol
col16: GeneID: NCBI GeneID, for those RefSeq genomes included in NCBI's Gene 
       resource
col17: locus_tag
col18: feature_interval_length: sum of the lengths of all intervals for the 
       feature (i.e. the length without introns for a joined feature)
col19: product_length: length of the product corresponding to the 
       accession.version in column 11. Protein product lengths are in amino acid
       units, and do not include the stop codon which is included in column 18.
       Additionally, product_length may differ from feature_interval_length if 
       the product contains sequence differences vs. the genome, as found for 
       some RefSeq transcript and protein products based on mRNA sequences and 
       also for INSDC proteins that are submitted to correct genome 
       discrepancies.
col20: attributes: semi-colon delimited list of a controlled set of qualifiers.
       The list currently includes:
       partial, pseudo, pseudogene, ribosomal_slippage, trans_splicing, 
       anticodon=NNN (for tRNAs), old_locus_tag=XXX 


*_genomic_gaps.txt.gz
---------------------
Tab-delimited text file reporting the coordinates of all gaps in the top-level 
genomic sequences. The gaps reported include gaps specified in the AGP files, 
gaps annotated on the component sequences, and any other run of 10 or more Ns 
in the sequences. Gap types are reported using the International Nucleotide 
Sequence Database Collaboration feature table terms with spaces replaced by 
underscores, see: http://www.insdc.org/files/feature_table.html

The file is tab delimited (including a #header) with the following columns:
col 1: sequence accession.version
col 2: gap start position (1-based)
col 3: gap stop position (1-based)
col 4: gap_length
col 5: gap_type. 
       One of: centromere, heterochromatin, short_arm, telomere,
       between_scaffolds, within_scaffold, repeat_between_scaffolds, 
       repeat_within_scaffold, contamination, unknown
col 6: linkage_evidence. 
       Gaps of type within_scaffold or repeat_within_scaffold have one or more 
       of the following types of linkage evidence: paired-ends, pcr, 
       proximity_ligation, align_genus, align_xgenus, align_trnscpt, 
       within_clone, clone_contig, map, strobe, unspecified, 
       inferred_from_sequence. 
       Multiple lines of linkage evidence are separated by a ';' delimiter. 
       Gaps of type contamination have unspecified as the linkage evidence.
       All other gap types have 'na' in the linkage evidence column. 
       See: http://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ 


annotation_hashes.txt
---------------------
Tab-delimited text file reporting hash values and change dates for specific 
details of the annotation. Hashes are computed based on the underlying data in 
ASN.1 format, and thus aren't affected by changes in file formats. In contrast,
the checksums reported in the md5checksums.txt file will change with any change
to the files, including file formats and differences in gzip compression. The
hashes are useful to monitor for when annotation has changed in a way that is 
significant for a particular use case and warrants downloading the updated 
records.

The file is tab delimited (including a #header) with the following columns:
col 1: Assembly accession: accession.version
col 2: Descriptors hash: hash of all descriptors on top-level sequence records,
       including BioSource, molinfo, user objects, publications, and dates
col 3: Descriptors last changed: date and time of the last change to any 
       descriptors
col 4: Features hash: hash of all features annotated on the assembly, including
       both locations and qualifiers stored directly on the genome records. For
       RefSeq genomes annotated with WP proteins and some other cases, protein
       product names aren't stored on the genome records and thus changes in 
       protein names do not alter the features hash.
col 5: Features last changed: date and time of the last change to any features
col 6: Locations hash: hash of just the locations of all features annotated on
       the assembly.
col 7: Locations last changed: date and time of the last change to any feature
       locations
col 8: Protein names hash: hash of the protein names for all CDS features 
       annotated on the assembly.
col 9: Protein names last changed: date and time of the last change to any 
       protein names.

Example use cases:
  A change in the Locations hash indicates that at least one feature has been 
     added, removed, or had its location altered.
  A change in the Features hash but not the Locations hash implies that only
     feature qualifiers have changed, such as names or db_xrefs.
  A change in the Protein names hash indicates that at least one protein name
     has changed compared to the previous files provided on the genomes FTP 
     site. Note for RefSeq prokaryotic genomes, protein names are updated 
     continuously but files on the FTP site are only refreshed intermittently
     to minimize churn.
  A change in the Descriptors hash but not the Features hash implies that only
     record metadata has been touched, such as the addition of a publication.

NOTE: currently the descriptors hash values are not stable due to a bug.


assembly_status.txt
------------------
A text file reporting the current status of the version of the assembly for 
which data is provided. Any assembly anomalies are also reported. Lines have the
format tag=value.

First line: status=<value> 
  where <value> is one of latest, replaced or suppressed
Second line (if any): assembly anomaly=<value>
  where value is a comma separated list of assembly anomalies as described in
  the "Anomalous assemblies" section of this web page:
  https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/


*_compare_prev.txt.gz
--------------------- 
The annotation produced for this release was compared to the annotation in the 
previous release. Scores for pairs of best-mapping current and previous gene and
transcript features were calculated based on overlap in exon sequence and 
matches in exon boundaries. Pairs of current and previous features were 
categorized based on these scores, whether they are reciprocal best matches, and
changes in attributes (gene biotype, completeness, etc.). If the assembly was 
updated between the two releases, alignments between the current and the 
previous assembly were used to match the current and previous gene and 
transcript features in aligning regions.

col 1: gene category: categorization of the difference between the gene in the 
       current annotation and the gene it maps best to in the previous 
       annotation: 
       * Changed locus ID: new feature identifier
       * Merged: the previous feature matches only part of the current feature
       * Split: the current feature matches only part of the previous feature
       * Changed locus type: difference in the gene biotype (coding vs. 
         non-coding, pseudogene vs. coding, etc...) - applies to genes only
       * Current-novel: Current feature has no matching previous feature
       * Previous-novel: previous feature has no matching current feature  
       * Current-other: current feature can't be matched unambiguously to 
         previous feature
       * Previous-other: previous feature can't be matched unambiguously to 
         current feature
       * Current-unmapped: current feature location can't be mapped to the 
         previous assembly
       * Previous-unmapped: previous feature location can't be mapped to the 
         current assembly
       * Moved: feature ID found on both current and previous but not placed 
         on regions aligned to each other by assembly-assembly alignment
       * Identical: identical exon boundaries in current and previous
       * Variant: alternative variant in current not in previous - applies 
         to transcripts only
       * Change in exception: Exceptions are added to RNA features when the 
         RefSeq transcript sequence doesn't match the conceptual sequence from
         the genome due to the presence of mismatches, indels, or additional 
         sequence, or in some other cases of unusual biology like ribosomal 
         slippage. This category reports when the current or previous RefSeq 
         transcript sequence was annotated with an exception and the matched 
         transcript does not. - applies to transcripts only
       * Similar: highly similar features, with support scores of 0.66 or more
         (on a scale of 0 to 1) on both sides of the comparison. The support 
         score is derived from a combination of matching exon boundaries and 
         sequence overlap. - applies to genes and non-coding transcripts only
       * Similar, change in CDS: support scores of 0.66 or more on both sides 
         AND the change affects the CDS - applies to coding transcripts only
       * Similar, change in UTR only: support scores of 0.66 or more on both 
         sides AND the change affects UTRs only (not the CDS) - applies to 
         coding transcripts only
       * Changed feature type: difference in the feature type - applies to 
         coding transcripts only
       * Changed completeness: feature is partial in current and complete in 
         the previous or vice versa - applies to genes only
       * Changed substantially: low similarity feature with support scores 
         below 0.66 on one or both sides. 
       * NA: No transcript feature associated with the current or previous gene 
         (pseudogene)
       * Previous-variant: alternative variant in previous not in current 
         - applies to transcripts only
       * Other: complex cases not fitting in categories above

       Notes:
       1. A gene may be categorized as not 'Identical' if the boundaries of a 
       gene are unchanged in the two annotation releases, but the children 
       exons have changed. This situation is frequent when a different set of 
       alternative variants are predicted for a gene. 
       2. Since a transcript feature on one side of the comparison may overlap
       multiple transcripts in the other side, a transcript may appear on 
       multiple lines in the report.
col 2: current GeneID: Gene database identifier for the current gene 
col 3: current gene biotype: attribute computed on gene features based on the 
       set of child features to indicate the overall biotype for the gene 
       annotation at this location. See list of possible biotypes in 
       ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt
col 4: current assembly unit: Assembly unit on which the current gene is 
       annotated   
col 5: current genomic accession: RefSeq accession.version of the sequence on 
       which the current gene is annotated
col 6: current gene range: start and stop coordinates of the current gene on the
       genomic sequence. 
col 7: orientation of the current gene on the genomic sequence
col 8-13: same description as for col 2-7 but for the previous gene that maps 
       best to the current gene.
col14: transcript category: categorization of the mapping between the current 
       and the previous annotated transcript. See description for column 1
col15: current transcript accession: RefSeq accession.version of the transcript
       annotated on the current assembly.
col16: current protein accession: RefSeq accession.version of the protein
       annotated on the current assembly.
col17: current transcript range: start and stop coordinates of the current 
       transcript on the genomic sequence.
col18-20: same description as for col 15-17 but for the previous transcript
       that maps best to the current transcript.
 
*_compare_prev.gbp.gz
---------------------
See the description for _compare_prev.txt.gz for a description of the content.
To visualize the differences in annotation in Genome WorkBench 
(https://www.ncbi.nlm.nih.gov/tools/gbench/):
Unzip *_compare_previous.gbp.gz
Load *_compare_previous.gbp to GenomeWorkBench 
Open the comparison report *_.compare_prev.table.asn as a Generic Table View
Use search and filter on terms in the table to find a specific gene, sequence 
or category
Right-click or double-click on a row and open the 'Graphical Sequence View' for
a graphical view of the gene on the current or the previous scaffold
Make sure 'Project features for aligned sequences' is enabled in the Alignments
track if the comparison is across two different versions of an assembly.

________________________________________________________________________________
National Center for Biotechnology Information (NCBI)
National Library of Medicine
National Institutes of Health
8600 Rockville Pike
Bethesda, MD 20894, USA
tel: (301) 496-2475
fax: (301) 480-9241
e-mail: info@ncbi.nlm.nih.gov
________________________________________________________________________________
