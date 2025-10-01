version 1.0

task partition_by_chromosome_task {
    input {
        File? inputBAM
        File? genome_fasta
        File? annot_gtf
        String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."

        String docker
    }

    command <<<
        set -ex
        set -o pipefail

        mkdir -p split_bams
        mkdir -p split_fastas
        mkdir -p split_gtfs
      
        if [ -f  "~{inputBAM}" ] && [ ! -f "~{inputBAM}.bai" ]; then
            samtools index ~{inputBAM}
        fi
        
        for chr in ~{chromosomes_want_partitioned}; do

            if [ -f "~{inputBAM}" ]; then
                samtools view -b ~{inputBAM} $chr > split_bams/$chr.bam
            fi

            if [ -f "~{genome_fasta}" ]; then
                samtools faidx ~{genome_fasta} $chr > split_fastas/$chr.genome.fasta
            fi
        
            if [ -f "~{annot_gtf}" ]; then
                cat ~{annot_gtf} | perl -lane 'if ($F[0] eq "'$chr'") { print; }' > split_gtfs/$chr.annot.gtf
            else        
                echo "# no gtf records" > split_gtfs/$chr.annot.gtf
            fi
        
        done
        
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_fastas/*.genome.fasta")
        Array[File] chromosomeGTFs = glob("split_gtfs/*.annot.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + ceil((size(inputBAM, "GB") + size(genome_fasta, "GB") + size(annot_gtf, "GB")) * 2.2 + 5) + " SSD"
    }
}


workflow partition_by_chromosome {
  input {
        File? inputBAM
        File? genome_fasta
        File? annot_gtf
        String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."

        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"

    }

    call partition_by_chromosome_task {
        input:
          inputBAM=inputBAM,
          genome_fasta=genome_fasta,
          annot_gtf=annot_gtf,
          chromosomes_want_partitioned=chromosomes_want_partitioned,
          docker = docker
      
    }

    output {
        Array[File] chromosomeBAMs = partition_by_chromosome_task.chromosomeBAMs
        Array[File] chromosomeFASTAs = partition_by_chromosome_task.chromosomeFASTAs
        Array[File] chromosomeGTFs = partition_by_chromosome_task.chromosomeGTFs
    }    
}


