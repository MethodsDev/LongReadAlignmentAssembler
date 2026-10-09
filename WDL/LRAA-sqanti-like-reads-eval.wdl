version 1.0

# Classifies long reads (a BAM) or isoforms (a GTF) into SQANTI-like categories.
#
# A BAM is classified in SHARDS, in three phases:
#
#   1. plan_shards    one small task: from the index alone, groups the small contigs
#                     in header order and cuts each contig above max_reads_per_shard
#                     into position ranges of similar read counts, plus one shard for
#                     the unplaced reads (util/misc/plan_bam_shards.py). A read belongs
#                     to the range holding its start, so none is classified twice, and
#                     no cut needs a gap in coverage: a read is classified on its own
#                     against its contig's annotation.
#   2. classify_shard one task per shard. Gets the whole BAM and its index, reads only
#                     its contigs or range from it (samtools view, as SAM into the
#                     classifier) and keeps only its contigs' lines of the annotation.
#                     The price of not cutting the BAM once is that every shard
#                     localizes all of it.
#   3. gather_shards  one task: appends the per-read tables and the tagged BAMs in
#                     shard order and sums the category counts, which reproduces what
#                     a single run over the whole BAM writes.
#
# A GTF of isoforms is small and is classified by one task.

workflow LRAA_sqanti_like_reads_eval_wf {
    input {
        String sample_id
        File ref_annot_GTF

        File? input_BAM
        File? input_BAI
        File? input_GTF

        # Named docker_sc, not docker: the R plotting step needs lraa-sc. Renamed 2026-10-06
        # so a saved Terra value for the old `docker` input (lraa:latest, the v0.17.7
        # default) cannot carry over and run this on lraa-core, which has no R.
        String docker_sc = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"

        # Contigs are grouped, in header order, until a group holds about this many
        # records; a contig above it is cut into ranges of about this many. On the
        # 46 million read test bam the largest contigs hold 3 to 4 million reads, so
        # they become 3 shards of 1.1 to 1.4 million (measured), and the many small
        # alt/unplaced contigs share a few. A shard is one task, so a very low value
        # only adds per-task overhead.
        Int max_reads_per_shard = 1500000

        # A shard reads its reads as a stream, so it classifies on one core whatever
        # its size; the cores only serve samtools and the compression of its output.
        # 8 GB is headroom, not a model: the single-task run of the whole 26 GB bam
        # peaked at 2.6 GiB with all 8 contig workers running at once (Terra
        # monitoring.log), so one worker needs a fraction of that.
        Int plan_cpu = 4
        Int shard_cpu = 4
        Int shard_memory_GB = 8

        Int preemptible_tries = 3

        # SSD by default: Terra's standard HDD throughput scales with disk size, so a
        # few hundred GB of HDD gives only tens of MB/s.
        String disk_type = "SSD"
    }

    if (defined(input_BAM)) {
        call plan_shards {
            input:
                input_BAM = select_first([input_BAM]),
                input_BAI = input_BAI,
                max_reads_per_shard = max_reads_per_shard,
                docker = docker_sc,
                cpu = plan_cpu,
                preemptible_tries = preemptible_tries,
                disk_type = disk_type
        }

        scatter (i in range(length(plan_shards.shard_specs))) {
            call classify_shard {
                input:
                    shard_name = "shard_" + i,
                    shard_spec = plan_shards.shard_specs[i],
                    input_BAM = select_first([input_BAM]),
                    input_BAI = plan_shards.bai,
                    ref_annot_GTF = ref_annot_GTF,
                    docker = docker_sc,
                    cpu = shard_cpu,
                    memory_GB = shard_memory_GB,
                    preemptible_tries = preemptible_tries,
                    disk_type = disk_type
            }
        }

        call gather_shards {
            input:
                sample_id = sample_id,
                shard_tsvs = classify_shard.iso_cats_tsv,
                shard_bams = classify_shard.iso_cats_bam,
                shard_summaries = classify_shard.summary_counts_tsv,
                docker = docker_sc,
                disk_type = disk_type
        }
    }

    # input_BAM wins when both are given.
    if (!defined(input_BAM) && defined(input_GTF)) {
        call classify_gtf {
            input:
                sample_id = sample_id,
                ref_annot_GTF = ref_annot_GTF,
                input_GTF = select_first([input_GTF]),
                docker = docker_sc,
                disk_type = disk_type
        }
    }

    output {
        File iso_cats_tsv = select_first([gather_shards.iso_cats_tsv, classify_gtf.iso_cats_tsv])
        File iso_cats_summary_counts_tsv = select_first([gather_shards.iso_cats_summary_counts_tsv, classify_gtf.iso_cats_summary_counts_tsv])
        File iso_cats_summary_counts_pdf = select_first([gather_shards.iso_cats_summary_counts_pdf, classify_gtf.iso_cats_summary_counts_pdf])
        File? iso_cats_bam = gather_shards.iso_cats_bam
    }
}


task plan_shards {
    input {
        File input_BAM
        File? input_BAI

        Int max_reads_per_shard
        String docker
        Int cpu = 4
        Int memory_GB = 8
        Int preemptible_tries
        String disk_type
    }

    Int disk_GB = ceil(1.2 * size(input_BAM, "GB") + 20)
    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu_tier = if cpu <= 4 then 4
        else if cpu <= 8 then 8
        else if cpu <= 16 then 16
        else if cpu <= 30 then 30
        else if cpu <= 60 then 60
        else if cpu <= 90 then 90
        else if cpu <= 180 then 180
        else 360
    Int c3d_mem_tier = if memory_GB <= 32 then 4
        else if memory_GB <= 64 then 8
        else if memory_GB <= 128 then 16
        else if memory_GB <= 240 then 30
        else if memory_GB <= 480 then 60
        else if memory_GB <= 720 then 90
        else if memory_GB <= 1440 then 180
        else 360
    Int c3d_effective_cpu = if c3d_cpu_tier >= c3d_mem_tier then c3d_cpu_tier else c3d_mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int c3d_highcpu_ram = if c3d_effective_cpu == 4 then 8
        else if c3d_effective_cpu == 8 then 16
        else if c3d_effective_cpu == 16 then 32
        else if c3d_effective_cpu == 30 then 59
        else if c3d_effective_cpu == 60 then 118
        else if c3d_effective_cpu == 90 then 177
        else if c3d_effective_cpu == 180 then 354
        else 708
    String c3d_machine_type = if memory_GB <= c3d_highcpu_ram
        then "c3d-highcpu-${c3d_effective_cpu}"
        else if memory_GB <= c3d_effective_cpu * 4
        then "c3d-standard-${c3d_effective_cpu}"
        else "c3d-highmem-${c3d_effective_cpu}"

    command <<<
        set -euo pipefail

        # The bam and bai can be localized to different directories; link both here
        # and index only if no bai was given.
        ln -s ~{input_BAM} input.bam
        if [[ "~{input_BAI}" != "" ]]; then
            ln -s ~{input_BAI} input.bam.bai
        else
            samtools index -@ ~{c3d_effective_cpu} input.bam
        fi

        # A copy, not the link: the index is an output, and a link to an input cannot be
        # delocalized.
        cp -L input.bam.bai plan.bam.bai

        # One line per shard, see the script for the format: whole contigs (small ones
        # grouped in header order), or a position range of one contig for a contig
        # above max_reads_per_shard, cut from the index alone. The unplaced reads, if
        # any, are the last shard. Shards stay in the bam's own order, which is what
        # lets the gather append them and get the order of a single run.
        "$(dirname "$(which SQANTI-like_cats_for_reads_or_isoforms.py)")/misc/plan_bam_shards.py" \
            --bam input.bam --bai input.bam.bai \
            --max_reads_per_shard ~{max_reads_per_shard} \
            > shard_specs.txt
    >>>

    output {
        Array[String] shard_specs = read_lines("shard_specs.txt")
        # the index the shards need, whether it was given or made here
        File bai = "plan.bam.bai"
    }

    runtime {
        docker: docker
        predefinedMachineType: c3d_machine_type
        bootDiskSizeGb: 50
        preemptible: preemptible_tries
        disks: "local-disk " + disk_GB + " " + disk_type
    }
}


task classify_shard {
    input {
        String shard_name
        # one line of the plan: "contigs<TAB>name name ..." (* is the unplaced reads) or
        # "range<TAB>name<TAB>start<TAB>end", the reads of one contig that start in it
        String shard_spec
        File input_BAM
        File input_BAI
        File ref_annot_GTF

        String docker
        Int cpu
        Int memory_GB
        Int preemptible_tries
        String disk_type
    }

    # Each shard localizes the whole bam and reads only its contigs from it, through
    # the index. The output is one shard's share: the tagged bam and the tables, both
    # at most the size of the bam.
    Int disk_GB = ceil(1.5 * size(input_BAM, "GB") + 3 * size(ref_annot_GTF, "GB") + 20)
    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu_tier = if cpu <= 4 then 4
        else if cpu <= 8 then 8
        else if cpu <= 16 then 16
        else if cpu <= 30 then 30
        else if cpu <= 60 then 60
        else if cpu <= 90 then 90
        else if cpu <= 180 then 180
        else 360
    Int c3d_mem_tier = if memory_GB <= 32 then 4
        else if memory_GB <= 64 then 8
        else if memory_GB <= 128 then 16
        else if memory_GB <= 240 then 30
        else if memory_GB <= 480 then 60
        else if memory_GB <= 720 then 90
        else if memory_GB <= 1440 then 180
        else 360
    Int c3d_effective_cpu = if c3d_cpu_tier >= c3d_mem_tier then c3d_cpu_tier else c3d_mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int c3d_highcpu_ram = if c3d_effective_cpu == 4 then 8
        else if c3d_effective_cpu == 8 then 16
        else if c3d_effective_cpu == 16 then 32
        else if c3d_effective_cpu == 30 then 59
        else if c3d_effective_cpu == 60 then 118
        else if c3d_effective_cpu == 90 then 177
        else if c3d_effective_cpu == 180 then 354
        else 708
    String c3d_machine_type = if memory_GB <= c3d_highcpu_ram
        then "c3d-highcpu-${c3d_effective_cpu}"
        else if memory_GB <= c3d_effective_cpu * 4
        then "c3d-standard-${c3d_effective_cpu}"
        else "c3d-highmem-${c3d_effective_cpu}"

    command <<<
        set -euo pipefail

        ln -s ~{input_BAM} input.bam
        ln -s ~{input_BAI} input.bam.bai

        IFS=$'\t' read -r kind field1 field2 field3 <<< '~{shard_spec}'

        # Each contig as {name} so samtools takes it literally, which contigs like
        # HLA-A*01:01 (GRCh38 alt haplotypes) need; * is the unplaced reads. Read into
        # an array rather than expanded unquoted, which would glob the *.
        regions=()
        min_pos=0
        if [ "$kind" = "range" ]; then
            names=("$field1")
            regions+=("{$field1}:$field2-$field3")
            # The query returns the reads overlapping the range, including ones that
            # began before it; those belong to the range before. A read belongs to the
            # range holding its start, so no read is classified twice.
            min_pos=$field2
        else
            read -r -a names <<< "$field1"
            for c in "${names[@]}"; do
                if [ "$c" = "*" ]; then regions+=("*"); else regions+=("{$c}"); fi
            done
        fi

        # The reference lines of this shard's contigs, so the script does not parse the
        # whole annotation. Empty for the unplaced reads, which are not classified.
        # (A range keeps its whole contig: a read is classified against all of it.)
        printf '%s\n' "${names[@]}" > shard_contigs.txt
        zcat -f ~{ref_annot_GTF} | awk -F'\t' '
            NR == FNR { want[$1] = 1; next }
            /^#/ { next }
            ($1 in want)
        ' shard_contigs.txt - > shard.gtf

        # The reads reach the script as SAM on stdin, header included, so no slice of
        # the bam is written. A stream has no index, so the contigs of the shard are
        # classified one after the other.
        samtools view -h input.bam "${regions[@]}" | \
            awk -F'\t' -v min_pos="$min_pos" '/^@/ || $4 >= min_pos' | \
            SQANTI-like_cats_for_reads_or_isoforms.py \
                --ref_gtf shard.gtf \
                --output_prefix ~{shard_name} \
                --input_bam - \
                --gzip_tsv --no_tsv_header --no_plot
    >>>

    output {
        File iso_cats_tsv = "~{shard_name}.iso_cats.tsv.gz"
        File iso_cats_bam = "~{shard_name}.iso_cats.bam"
        File summary_counts_tsv = "~{shard_name}.iso_cats.summary_counts.tsv"
    }

    runtime {
        docker: docker
        predefinedMachineType: c3d_machine_type
        bootDiskSizeGb: 50
        preemptible: preemptible_tries
        disks: "local-disk " + disk_GB + " " + disk_type
    }
}


task gather_shards {
    input {
        String sample_id
        Array[File] shard_tsvs
        Array[File] shard_bams
        Array[File] shard_summaries

        String docker
        Int cpu = 4
        Int memory_GB = 8
        String disk_type
        # Not preemptible: one task holding every shard's output, so a preemption
        # discards the whole gather.
        Int preemptible_tries = 0
    }

    # the parts, and the same bytes again in the appended table and bam
    Int disk_GB = ceil(2.2 * (size(shard_tsvs, "GB") + size(shard_bams, "GB")) + 20)
    
    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu_tier = if cpu <= 4 then 4
        else if cpu <= 8 then 8
        else if cpu <= 16 then 16
        else if cpu <= 30 then 30
        else if cpu <= 60 then 60
        else if cpu <= 90 then 90
        else if cpu <= 180 then 180
        else 360
    Int c3d_mem_tier = if memory_GB <= 32 then 4
        else if memory_GB <= 64 then 8
        else if memory_GB <= 128 then 16
        else if memory_GB <= 240 then 30
        else if memory_GB <= 480 then 60
        else if memory_GB <= 720 then 90
        else if memory_GB <= 1440 then 180
        else 360
    Int c3d_effective_cpu = if c3d_cpu_tier >= c3d_mem_tier then c3d_cpu_tier else c3d_mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int c3d_highcpu_ram = if c3d_effective_cpu == 4 then 8
        else if c3d_effective_cpu == 8 then 16
        else if c3d_effective_cpu == 16 then 32
        else if c3d_effective_cpu == 30 then 59
        else if c3d_effective_cpu == 60 then 118
        else if c3d_effective_cpu == 90 then 177
        else if c3d_effective_cpu == 180 then 354
        else 708
    String c3d_machine_type = if memory_GB <= c3d_highcpu_ram
        then "c3d-highcpu-${c3d_effective_cpu}"
        else if memory_GB <= c3d_effective_cpu * 4
        then "c3d-standard-${c3d_effective_cpu}"
        else "c3d-highmem-${c3d_effective_cpu}"

    command <<<
        set -euo pipefail

        # The shards are given in shard order, which is contig header order. The two
        # appends are independent copies (no recompression), so they run side by side.
        # `wait <pid>` returns that job's status, so a failed copy fails the task.
        #
        # Table: header, then each shard's table as it is. Appended gzip members are one
        # valid gzip stream (zcat, Python and R read it as one).
        (
            python3 -c '
import gzip, importlib.machinery, shutil, sys
sq = importlib.machinery.SourceFileLoader("sq", shutil.which("SQANTI-like_cats_for_reads_or_isoforms.py")).load_module()
with gzip.open("~{sample_id}.iso_cats.tsv.gz", "wt", compresslevel=6) as ofh:
    ofh.write("\t".join(sq.BAM_TSV_FIELDNAMES) + "\n")
'
            while read -r f; do
                cat "$f" >> ~{sample_id}.iso_cats.tsv.gz
            done < ~{write_lines(shard_tsvs)}
        ) &
        tsv_pid=$!

        # Bam: the shards cover disjoint contigs, in header order, so appending them
        # keeps the bam coordinate sorted; samtools cat copies the compressed blocks.
        samtools cat --no-PG -b ~{write_lines(shard_bams)} -o ~{sample_id}.iso_cats.bam &
        bam_pid=$!

        wait $tsv_pid
        wait $bam_pid

        # Summed in slice order, so categories stay in the order a single run first sees them.
        awk -F'\t' -v OFS='\t' '
            FNR == 1 { next }
            !($1 in count) { order[++n] = $1 }
            { count[$1] += $2 }
            END {
                print "Category", "Count"
                for (i = 1; i <= n; i++) print order[i], count[order[i]]
            }
        ' $(cat ~{write_lines(shard_summaries)}) > ~{sample_id}.iso_cats.summary_counts.tsv

        "$(dirname "$(which SQANTI-like_cats_for_reads_or_isoforms.py)")/misc/plot_SQANTI_cats.Rscript" \
            ~{sample_id}.iso_cats.summary_counts.tsv ~{sample_id}.iso_cats.summary_counts.pdf

        gzip ~{sample_id}.iso_cats.summary_counts.tsv
    >>>

    output {
        File iso_cats_tsv = "~{sample_id}.iso_cats.tsv.gz"
        File iso_cats_summary_counts_tsv = "~{sample_id}.iso_cats.summary_counts.tsv.gz"
        File iso_cats_summary_counts_pdf = "~{sample_id}.iso_cats.summary_counts.pdf"
        File iso_cats_bam = "~{sample_id}.iso_cats.bam"
    }

    runtime {
        docker: docker
        predefinedMachineType: c3d_machine_type
        bootDiskSizeGb: 50
        preemptible: preemptible_tries
        disks: "local-disk " + disk_GB + " " + disk_type
    }
}


task classify_gtf {
    input {
        String sample_id
        File ref_annot_GTF
        File input_GTF

        String docker
        Int cpu = 4
        Int memory_GB = 16
        String disk_type
    }

    Int disk_GB = ceil(4 * size(input_GTF, "GB") + 4 * size(ref_annot_GTF, "GB") + 50)
    
    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu_tier = if cpu <= 4 then 4
        else if cpu <= 8 then 8
        else if cpu <= 16 then 16
        else if cpu <= 30 then 30
        else if cpu <= 60 then 60
        else if cpu <= 90 then 90
        else if cpu <= 180 then 180
        else 360
    Int c3d_mem_tier = if memory_GB <= 32 then 4
        else if memory_GB <= 64 then 8
        else if memory_GB <= 128 then 16
        else if memory_GB <= 240 then 30
        else if memory_GB <= 480 then 60
        else if memory_GB <= 720 then 90
        else if memory_GB <= 1440 then 180
        else 360
    Int c3d_effective_cpu = if c3d_cpu_tier >= c3d_mem_tier then c3d_cpu_tier else c3d_mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int c3d_highcpu_ram = if c3d_effective_cpu == 4 then 8
        else if c3d_effective_cpu == 8 then 16
        else if c3d_effective_cpu == 16 then 32
        else if c3d_effective_cpu == 30 then 59
        else if c3d_effective_cpu == 60 then 118
        else if c3d_effective_cpu == 90 then 177
        else if c3d_effective_cpu == 180 then 354
        else 708
    String c3d_machine_type = if memory_GB <= c3d_highcpu_ram
        then "c3d-highcpu-${c3d_effective_cpu}"
        else if memory_GB <= c3d_effective_cpu * 4
        then "c3d-standard-${c3d_effective_cpu}"
        else "c3d-highmem-${c3d_effective_cpu}"

    command <<<
        set -ex

        SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf ~{ref_annot_GTF} --output_prefix ~{sample_id} --input_gtf ~{input_GTF} --gzip_tsv

        # the per-feature table is already gzipped (--gzip_tsv); only the summary is left
        gzip ~{sample_id}.iso_cats.summary_counts.tsv
    >>>

    output {
        File iso_cats_tsv = "~{sample_id}.iso_cats.tsv.gz"
        File iso_cats_summary_counts_tsv = "~{sample_id}.iso_cats.summary_counts.tsv.gz"
        File iso_cats_summary_counts_pdf = "~{sample_id}.iso_cats.summary_counts.pdf"
    }

    runtime {
        docker: docker
        predefinedMachineType: c3d_machine_type
        bootDiskSizeGb: 50
        disks: "local-disk " + disk_GB + " " + disk_type
    }
}
