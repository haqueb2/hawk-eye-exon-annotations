configfile: "config.yaml"

rule all:
    input:
        "1-Parse ExonCalculator/HAWKEYE_Database_Parse_ExonCalculator.csv",
        "2-ClinVar/clinvar_LPP_variants.csv",
        "2-ClinVar/clinvar_LPP_variants_in_HAWKEYE_exons.csv",
        "2-ClinVar/HAWKEYE_Database_ClinVar_Counts.csv",
        "3-Premature Stop/HAWKEYE_Database_Exon_Features.csv",
        "4-PhyloP/HAWKEYE_Database_PhyloP_Exon_Conservation.csv",
        "5-UniProt/HAWKEYE_Database_UniProt_All_Protein_Domain_Features.csv",
        "6-InterPro/HAWKEYE_Database_InterPro_Features.csv",
        "7-GTEx/HAWKEYE_Database_GTEx_Expression.csv",
        "8-PSI/HAWKEYE_Database_PSI.csv",
        "9-Domain Counting/HAWKEYE_Database_Domain_Repeat_Counts.csv"


rule parse_exoncalculator:
    input:
        exon_length="1-Parse ExonCalculator/ExonLength_output.csv",
        cds="1-Parse ExonCalculator/CDS_output.csv"
    output:
        "1-Parse ExonCalculator/HAWKEYE_Database_Parse_ExonCalculator.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '1-Parse ExonCalculator/parse_exoncalculator.R' 'input_file=\"{input.exon_length}\"' 'cds_file=\"{input.cds}\"' 'output_file=\"{output}\"'"


rule clinvar_annotation:
    input:
        exon="1-Parse ExonCalculator/HAWKEYE_Database_Parse_ExonCalculator.csv",
        clinvar_summary="2-ClinVar/clinvar_db/variant_summary.txt",
        clinvar_inframe="2-ClinVar/clinvar_db/CLINVAR_20251106_INS_DEL_DUP_INV.csv"
    output:
        all_lpp="2-ClinVar/clinvar_LPP_variants.csv",
        lpp_exonic="2-ClinVar/clinvar_LPP_variants_in_HAWKEYE_exons.csv",
        exon_counts="2-ClinVar/HAWKEYE_Database_ClinVar_Counts.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '2-ClinVar/lpp_clinvar_annotation.R' 'exon_file=\"{input.exon}\"' 'clinvar_summary_file=\"{input.clinvar_summary}\"' 'clinvar_inframe_file=\"{input.clinvar_inframe}\"' 'output_all_lpp=\"{output.all_lpp}\"' 'output_lpp_exonic=\"{output.lpp_exonic}\"' 'output_exon_counts=\"{output.exon_counts}\"'"


rule exon_features:
    input:
        exon_counts="2-ClinVar/HAWKEYE_Database_ClinVar_Counts.csv",
        ccds_seq="3-Premature Stop/premature_stop/CCDS_nucleotide.current.fna",
        ccds_map="3-Premature Stop/premature_stop/CCDS2Sequence.current.txt",
        ccds_summary="3-Premature Stop/premature_stop/CCDS.current.txt"
    output:
        "3-Premature Stop/HAWKEYE_Database_Exon_Features.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '3-Premature Stop/exon_features.R' 'exon_file=\"{input.exon_counts}\"' 'ccds_sequence_file=\"{input.ccds_seq}\"' 'ccds_to_sequence_file=\"{input.ccds_map}\"' 'ccds_summary_file=\"{input.ccds_summary}\"' 'output_file=\"{output}\"'"


rule phylop_conservation:
    input:
        exon_features="3-Premature Stop/HAWKEYE_Database_Exon_Features.csv",
        bigwig="4-PhyloP/phylop/hg38.phyloP100way.bw"
    output:
        "4-PhyloP/HAWKEYE_Database_PhyloP_Exon_Conservation.csv"
    conda: "envs/py.yaml"
    shell:
        "python3 4-PhyloP/phylop_mean_conservation_score.py --input '{input.exon_features}' --bigwig '{input.bigwig}' --output '{output}'"


rule uniprot_features_step1:
    input:
        exon_features="4-PhyloP/HAWKEYE_Database_PhyloP_Exon_Conservation.csv"
    output:
        "5-UniProt/HAWKEYE_Database_UniProt_Features.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '5-UniProt/uniprot_features.R' 'input_file=\"{input.exon_features}\"' 'idmapping_file=\"5-UniProt/uniprot_files/HUMAN_9606_idmapping.dat\"' 'output_file=\"{output}\"'"


rule protein_domain_annotation:
    input:
        uniprot_features="5-UniProt/HAWKEYE_Database_UniProt_Features.csv"
    output:
        "5-UniProt/HAWKEYE_Database_UniProt_All_Protein_Domain_Features.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '5-UniProt/protein_domain_annotation.R' 'input_file=\"{input.uniprot_features}\"' 'output_file=\"{output}\"'"


rule uniprot_features_step2:
    input:
        exon_features="4-PhyloP/HAWKEYE_Database_PhyloP_Exon_Conservation.csv"
    output:
        "5-UniProt/HAWKEYE_Database_UniProt_Features_after_domains.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '5-UniProt/uniprot_features.R' 'input_file=\"{input.exon_features}\"' 'idmapping_file=\"5-UniProt/uniprot_files/HUMAN_9606_idmapping.dat\"' 'output_file=\"{output}\"'"


rule interpro_annotation:
    input:
        uniprot_all="5-UniProt/HAWKEYE_Database_UniProt_All_Protein_Domain_Features.csv"
    output:
        "6-InterPro/HAWKEYE_Database_InterPro_Features.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '6-InterPro/interpro_protein_domain_annotation.R' 'input_file=\"{input.uniprot_all}\"' 'output_file=\"{output}\"'"


rule gtex_expression:
    input:
        interpro="6-InterPro/HAWKEYE_Database_InterPro_Features.csv"
    output:
        "7-GTEx/HAWKEYE_Database_GTEx_Expression.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '7-GTEx/add_gtex_expression.R' 'input_file=\"{input.interpro}\"' 'output_file=\"{output}\"'"


rule psi:
    input:
        gtex="7-GTEx/HAWKEYE_Database_GTEx_Expression.csv"
    output:
        "8-PSI/HAWKEYE_Database_PSI.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '8-PSI/ensembl_psi.R' 'input_file=\"{input.gtex}\"' 'output_file=\"{output}\"'"


rule domain_counting:
    input:
        psi="8-PSI/HAWKEYE_Database_PSI.csv",
        uniprot_all="5-UniProt/HAWKEYE_Database_UniProt_All_Protein_Domain_Features.csv"
    output:
        "9-Domain Counting/HAWKEYE_Database_Domain_Repeat_Counts.csv"
    conda: "envs/r.yaml"
    shell:
        "Rscript '9-Domain Counting/uniprot_protein_domain_repeat_counting.R' 'input_file=\"{input.psi}\"' 'output_file=\"{output}\"' 'uniprot_all=\"{input.uniprot_all}\"'"
