#!/bin/bash

source /conda.sh
conda activate /conda_env

BAM_DIR="/rnaseq/aligned"
GTF_FILE="/UBA1.gtf"

OUTPUT_DIR="/uba1_exon_enrichment_fpkm"
mkdir -p "$OUTPUT_DIR"
OUTPUT_CSV="$OUTPUT_DIR/exons_fpkm.csv"

echo "Extracting exon coordinates from GTF..."
EXON_COORDS=$(awk '$3 == "exon" {print $1, $4, $5, $10}' "$GTF_FILE")

#Epected file pairs
declare -A FILE_PAIRS=(
    ["VEXAS_CD14"]="$BAM_DIR/VEXAS_CD14Aligned.sortedByCoord.out.bam $BAM_DIR/VEXAS_CD14SJ.out.tab"
    ["VEXAS_CD3"]="$BAM_DIR/VEXAS_CD3Aligned.sortedByCoord.out.bam $BAM_DIR/VEXAS_CD3SJ.out.tab"
    ["CTRL_CD14"]="$BAM_DIR/CTRL_CD14Aligned.sortedByCoord.out.bam $BAM_DIR/CTRL_CD14SJ.out.tab"
    ["CTRL_CD3"]="$BAM_DIR/CTRL_CD3Aligned.sortedByCoord.out.bam $BAM_DIR/CTRL_CD3SJ.out.tab"
)

# Ensure BAM files are indexed
for SAMPLE in "${!FILE_PAIRS[@]}"; do
    BAM_FILE="${FILE_PAIRS[$SAMPLE]%% *}"
    if [[ ! -f "$BAM_FILE.bai" ]]; then
        echo "⚠️ Index file for $BAM_FILE is missing. Exiting."
        exit 1
    fi
done

#Output header
echo "Sample,Order,Chromosome,Start,End,Strand,Motif,Annotated,Unique_Reads,Multi_Reads,Max_Splice_Score,Exon_Coords_Start,Exon_Coords_End,Reads,FPKM,Isoform_Exon_Overlap,Predicted_Exon" > "$OUTPUT_CSV"

#Process each sample explicitly
for SAMPLE in "${!FILE_PAIRS[@]}"; do
    BAM_FILE="${FILE_PAIRS[$SAMPLE]%% *}"
    SJ_FILE="${FILE_PAIRS[$SAMPLE]##* }"

    echo "----------------------------------"
    echo "Checking files for: $SAMPLE"
    echo "BAM File: $BAM_FILE"
    echo "SJ File: $SJ_FILE"

    if [[ ! -f "$BAM_FILE" || ! -f "$SJ_FILE" ]]; then
        echo "Missing BAM or SJ.out.tab for $SAMPLE. Skipping..."
        continue
    else
        echo "Both BAM and SJ.out.tab found!"
    fi
    echo "----------------------------------"

    #Process the SJ.out.tab file
    awk -v OFS="," -v exons="$EXON_COORDS" -v bam="$BAM_FILE" '
    BEGIN {
        split(exons, exon_list, "\n");
        for (i in exon_list) {
            split(exon_list[i], e, " ");
            exon_start[e[2]"-"e[3]] = e[1];
            exon_end[e[2]"-"e[3]] = e[1];
            exon_isoform[e[2]"-"e[3]] = e[4];
        }
        prev_junction_end = "";
    }
    ($1 == "chrX" && $2 >= 47190861 && $3 <= 47215128) {
        order = NR;
        junction_start = $2;
        junction_end = $3;
        
        #Check if a predicted exon exists between previous junction end and current start (within 3 bp)
        if (prev_junction_end != "" && prev_junction_end < junction_start - 3) {
            exon_coords_start = prev_junction_end;
            exon_coords_end = junction_start;
            
            #Find isoform overlap within 3 bp
            isoform_exon_overlap = "";
            for (i in exon_list) {
                split(exon_list[i], e, " ");
                if (e[2] >= exon_coords_start - 3 && e[3] <= exon_coords_end + 3) {  
                    isoform_exon_overlap = (isoform_exon_overlap == "") ? e[4] : isoform_exon_overlap "," e[4];
                }
            }

            #Count reads in exon region
            cmd = "samtools view -c -F 4 " bam " " exon_coords_start "-" exon_coords_end;
            cmd | getline read_count;
            close(cmd);

            #Get total reads mapped to UBA1 gene
            cmd2 = "samtools view -c -F 4 " bam " chrX:47190861-47215128";
            cmd2 | getline total_uba1_reads;
            close(cmd2);

            #Calculate exon length
            exon_length = exon_coords_end - exon_coords_start;

            #Compute FPKM normalization
            if (total_uba1_reads > 0 && exon_length > 0 && read_count > 0) {
                fpkm = (read_count * 1000000000) / (total_uba1_reads * exon_length);
            } else {
                fpkm = "N/A";
            }

            #Print predicted exon row
            print SAMPLE, order, $1, exon_coords_start, exon_coords_end, $4, $5, $6, "N/A", "N/A", "N/A", exon_coords_start, exon_coords_end, read_count, fpkm, isoform_exon_overlap, "Yes" >> OUTPUT_CSV;
        }
        
        #Update previous junction end
        prev_junction_end = junction_end;

        #Print the normal splice junction row
        print SAMPLE, order, $1, $2, $3, $4, $5, $6, $7, $8, $9, "N/A", "N/A", "N/A", "N/A", "N/A", "No" >> OUTPUT_CSV;
    }
    ' "$SJ_FILE"

    echo "Processed $SAMPLE and combined results in $OUTPUT_CSV"
done

echo "Annotation and read count analysis complete! Results saved in $OUTPUT_CSV"

