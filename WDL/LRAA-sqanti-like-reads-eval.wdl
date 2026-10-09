version 1.0

# Classifies long reads (a BAM) or isoforms (a GTF) into SQANTI-like categories.
#
# A BAM is classified in SHARDS, in three phases:
#
#   1. plan_shards    one small task: from the index alone, groups the small contigs
#                     in header order and cuts each contig above max_reads_per_shard
#                     into position ranges of similar read counts, plus one shard for
#                     the unplaced reads (util/misc/plan_bam_shards.py). It also cuts
#                     the BAM and the annotation into one file per shard. A read belongs
#                     to the range holding its start, so none is classified twice, and
#                     no cut needs a gap in coverage: a read is classified on its own
#                     against its contig's annotation.
#   2. classify_shard one task per shard. Gets only its own slice of the BAM (the plan
#                     task cut one per shard, so a shard localizes its share, not the
#                     whole BAM, which matters more the larger the BAM), indexes it, and
#                     runs one chunk per core: reads the chunk's contigs or range from
#                     the slice (samtools view, as SAM into the classifier).
#   3. gather_shards  one task: appends the per-read tables and the tagged BAMs of all
#                     chunks in order and sums the category counts, which reproduces
#                     what a single run over the whole BAM writes. The shards hand over
#                     their chunks' files as they are; appending them anywhere else
#                     would only be repeated here.
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

        # A CHUNK is what one process classifies: contigs are grouped, in header order,
        # until a group holds about this many records, and a contig above it is cut into
        # ranges of about this many (from the index alone). A SHARD is what one VM runs:
        # chunks_per_shard consecutive chunks, run in parallel on the VM's cores. Use
        # about twice as many chunks as cores: a shard starts its largest chunks first
        # and the smaller ones fill the cores as they free up, so no core waits at the
        # end on one big chunk. Measured (4 cores, 4 chunks of 1.1 to 1.3 million reads):
        # one chunk alone on a VM keeps 1.3 of 4 cores busy, four together keep 3.25 and
        # take 153 s of VM time against 390 s for four VMs; 8 chunks of about 600
        # thousand took 133 s.
        Int max_reads_per_chunk = 600000
        Int chunks_per_shard = 8

        # Each chunk is classified on one core, so a shard uses as many cores as it has
        # chunks in flight; the rest compress the tagged bam of a shard with few chunks.
        # 8 GB is headroom, not a model: the single-task run of the whole 26 GB bam
        # peaked at 2.6 GiB with all 8 contig workers running at once (Terra
        # monitoring.log), so one worker needs a fraction of that.
        # The plan task cuts one slice of the BAM per shard, one per core at a time, so
        # give it cores in proportion to the BAM: it decompresses and recompresses all
        # of it once.
        Int plan_cpu = 16
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
                ref_annot_GTF = ref_annot_GTF,
                max_reads_per_chunk = max_reads_per_chunk,
                max_reads_per_shard = max_reads_per_chunk * chunks_per_shard,
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
                    shard_bam = plan_shards.shard_bams[i],
                    shard_gtf = plan_shards.shard_gtfs[i],
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
                chunk_tsvs = flatten(classify_shard.chunk_tsvs),
                chunk_bams = flatten(classify_shard.chunk_bams),
                chunk_summaries = flatten(classify_shard.chunk_summaries),
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
        File ref_annot_GTF

        Int max_reads_per_chunk
        Int max_reads_per_shard
        String docker
        Int cpu = 4
        Int memory_GB = 8
        Int preemptible_tries
        String disk_type
    }

    # the bam, and its slices (level 1 compression, a little larger than the bam)
    Int disk_GB = ceil(2.5 * size(input_BAM, "GB") + 2.2 * size(ref_annot_GTF, "GB") + 20)
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

        # One line per shard, see the script for the format: its chunks, each whole
        # contigs (small ones grouped in header order) or a position range of one contig
        # above max_reads_per_chunk, cut from the index alone. The unplaced reads, if
        # any, are the last shard. Shards stay in the bam's own order, which is what
        # lets the gather append them and get the order of a single run.
        "$(dirname "$(which SQANTI-like_cats_for_reads_or_isoforms.py)")/misc/plan_bam_shards.py" \
            --bam input.bam --bai input.bam.bai --threads ~{c3d_effective_cpu} \
            --max_reads_per_chunk ~{max_reads_per_chunk} \
            --max_reads_per_shard ~{max_reads_per_shard} \
            > shard_specs.txt

        # One slice of the bam per shard: the union of its chunks' regions, taken through
        # the index, so the bam is read once in total and each shard localizes only its
        # share. Adjacent ranges of one contig are one query (a read spanning the
        # boundary is then in the slice once); each chunk still filters by start
        # position, so the slice needs nothing more than the reads overlapping its
        # regions. Contigs are written as {name} so samtools takes them literally, which
        # contigs like HLA-A*01:01 (GRCh38 alt haplotypes) need; * is the unplaced reads.
        mkdir shard_bams
        awk -F'\t' '{
            n = split($0, chunks, ",")
            regions = ""; prev = ""; run_start = ""; run_end = ""
            for (i = 1; i <= n; i++) {
                m = split(chunks[i], f, "\t")
                if (f[1] == "range") {
                    if (f[2] == prev && f[3] == run_end + 1) { run_end = f[4] }
                    else {
                        if (prev != "") regions = regions " \047{" prev "}:" run_start "-" run_end "\047"
                        prev = f[2]; run_start = f[3]; run_end = f[4]
                    }
                } else {
                    if (prev != "") regions = regions " \047{" prev "}:" run_start "-" run_end "\047"
                    prev = ""
                    k = split(f[2], names, " ")
                    for (j = 1; j <= k; j++) regions = regions " \047" (names[j] == "*" ? "*" : "{" names[j] "}") "\047"
                }
            }
            if (prev != "") regions = regions " \047{" prev "}:" run_start "-" run_end "\047"
            printf "samtools view -O bam,level=1 --no-PG -o shard_bams/shard_%05d.bam input.bam%s\n", NR - 1, regions
        }' shard_specs.txt > slice_commands.txt
        # independent reads of one bam through its index, as many at once as cores
        xargs -P ~{c3d_effective_cpu} -d '\n' -I CMD bash -c CMD < slice_commands.txt

        # One annotation file per shard, in one pass over the reference: the lines of the
        # contigs the shard's chunks need. A contig cut into ranges can be needed by two
        # shards, so a line goes to every shard that needs its contig. Every shard gets a
        # file, empty for the unplaced reads.
        mkdir shard_gtfs
        awk -F'\t' '{
            n = split($0, chunks, ",")
            for (i = 1; i <= n; i++) {
                m = split(chunks[i], f, "\t")
                if (f[1] == "range") { print f[2] "\t" NR - 1 }
                else { k = split(f[2], names, " "); for (j = 1; j <= k; j++) print names[j] "\t" NR - 1 }
            }
        }' shard_specs.txt | sort -u > contig_to_shards.tsv
        for i in $(seq 0 $(( $(wc -l < shard_specs.txt) - 1 ))); do
            : > "$(printf 'shard_gtfs/shard_%05d.gtf' "$i")"
        done
        zcat -f ~{ref_annot_GTF} | awk -F'\t' '
            NR == FNR { shards[$1] = shards[$1] " " $2; next }
            /^#/ { next }
            ($1 in shards) {
                k = split(shards[$1], ids, " ")
                for (i = 1; i <= k; i++) print > sprintf("shard_gtfs/shard_%05d.gtf", ids[i])
            }
        ' contig_to_shards.tsv -
    >>>

    output {
        Array[String] shard_specs = read_lines("shard_specs.txt")
        # index i is shard i's slice of the bam, and its annotation
        Array[File] shard_bams = glob("shard_bams/*.bam")
        Array[File] shard_gtfs = glob("shard_gtfs/*.gtf")
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
        # one line of the plan: chunks joined with commas, each
        # "contigs<TAB>name name ...<TAB>reads" (* is the unplaced reads) or
        # "range<TAB>name<TAB>start<TAB>end<TAB>reads", the reads of one contig that
        # start in it; reads is the estimated record count
        String shard_spec
        # this shard's slice of the bam: the reads overlapping its chunks' regions
        File shard_bam
        # the annotation lines of the contigs of this shard's chunks
        File shard_gtf

        String docker
        Int cpu
        Int memory_GB
        Int preemptible_tries
        String disk_type
    }

    # The slice, its index, and the chunks' outputs (the tagged bam and the tables, at
    # most the size of the slice).
    Int disk_GB = ceil(3 * size(shard_bam, "GB") + 3 * size(shard_gtf, "GB") + 20)
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

        ln -s ~{shard_bam} input.bam
        samtools index -@ ~{c3d_effective_cpu} input.bam

        # The chunks of this shard, in order: one spec file each, and the estimated read
        # count of each (the last field) for starting the largest first.
        mkdir chunks
        printf '%s' '~{shard_spec}' | tr ',' '\n' | awk -F'\t' '{
            print > sprintf("chunks/%05d.spec", NR - 1)
            print $NF "\t" NR - 1 > "chunk_sizes.tsv"
        }'
        num_chunks=$(ls chunks/*.spec | wc -l)

        # The shard's annotation, split once more: one file per contig any chunk needs, so
        # a chunk loads only its own contigs' lines (this file is small, the plan task
        # already cut it down to this shard's contigs).
        mkdir gtf
        cat chunks/*.spec | awk -F'\t' '$1 == "range" { print $2; next } { n = split($2, a, " "); for (i = 1; i <= n; i++) print a[i] }' | sort -u > contigs_needed.txt
        zcat -f ~{shard_gtf} | awk -F'\t' '
            NR == FNR { id[$1] = ++n; next }
            /^#/ { next }
            ($1 in id) { print > sprintf("gtf/%d.gtf", id[$1]) }
        ' contigs_needed.txt -
        # a contig with no annotation (and the unplaced reads) still gets a file
        n=0; while read -r c; do n=$((n + 1)); touch "gtf/$n.gtf"; done < contigs_needed.txt

        # One chunk: its reads as SAM, header included, straight into the classifier, so
        # no slice of the bam is written. Several of these run at once; each is one
        # process on one core, so the cores that no chunk needs compress the tagged bam.
        classify_chunk() {
            local j="$1"
            local prefix
            prefix=$(printf 'chunks/%05d' "$j")
            local out
            out=$(printf 'chunks/~{shard_name}.c%05d' "$j")
            local kind field1 field2 field3 estimate
            # contigs: kind, names, reads.  range: kind, name, start, end, reads.
            IFS=$'\t' read -r kind field1 field2 field3 estimate < "$prefix.spec"

            # Each contig as {name} so samtools takes it literally, which contigs like
            # HLA-A*01:01 (GRCh38 alt haplotypes) need; * is the unplaced reads. Read
            # into an array rather than expanded unquoted, which would glob the *.
            local regions=() names=() min_pos=0 c
            if [ "$kind" = "range" ]; then
                names=("$field1")
                regions+=("{$field1}:$field2-$field3")
                # The query returns the reads overlapping the range, including ones that
                # began before it; those belong to the range before. A read belongs to
                # the range holding its start, so no read is classified twice.
                min_pos=$field2
            else
                read -r -a names <<< "$field1"
                for c in "${names[@]}"; do
                    if [ "$c" = "*" ]; then regions+=("*"); else regions+=("{$c}"); fi
                done
            fi

            # this chunk's annotation: its contigs' files (a range keeps its whole
            # contig, since a read is classified against all of it)
            : > "$prefix.gtf"
            for c in "${names[@]}"; do
                cat "gtf/$(grep -n -x -F -- "$c" contigs_needed.txt | head -1 | cut -d: -f1).gtf" >> "$prefix.gtf"
            done

            samtools view -h input.bam "${regions[@]}" | \
                awk -F'\t' -v min_pos="$min_pos" '/^@/ || $4 >= min_pos' | \
                SQANTI-like_cats_for_reads_or_isoforms.py \
                    --ref_gtf "$prefix.gtf" \
                    --output_prefix "$out" \
                    --input_bam - \
                    --gzip_tsv --no_tsv_header --no_plot \
                    --bam_write_threads ${extra_threads}
        }
        export -f classify_chunk

        # Cores no chunk can use go to compressing the bam: 3 extra for a lone chunk on
        # 4 cores, none when there are as many chunks as cores.
        extra_threads=$(( ~{c3d_effective_cpu} / num_chunks - 1 ))
        if [ "$extra_threads" -lt 0 ]; then extra_threads=0; fi
        export extra_threads

        # As many chunks at once as cores, the largest first (by the plan's estimate) so
        # the smaller ones fill the cores as they free up. The order of execution does not
        # matter: the files carry the chunk number and the gather appends them in order.
        # A failed chunk fails the task.
        sort -k1,1nr chunk_sizes.tsv | cut -f2 \
            | xargs -P ~{c3d_effective_cpu} -I{} bash -c 'classify_chunk {}'
    >>>

    output {
        # one file per chunk, in chunk order: the file names carry the chunk number
        Array[File] chunk_tsvs = glob("chunks/*.iso_cats.tsv.gz")
        Array[File] chunk_bams = glob("chunks/*.iso_cats.bam")
        Array[File] chunk_summaries = glob("chunks/*.iso_cats.summary_counts.tsv")
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
        # every chunk's files, in the bam's order
        Array[File] chunk_tsvs
        Array[File] chunk_bams
        Array[File] chunk_summaries

        String docker
        Int cpu = 4
        Int memory_GB = 8
        String disk_type
        # Not preemptible: one task holding every shard's output, so a preemption
        # discards the whole gather.
        Int preemptible_tries = 0
    }

    # the parts, and the same bytes again in the appended table and bam
    Int disk_GB = ceil(2.2 * (size(chunk_tsvs, "GB") + size(chunk_bams, "GB")) + 20)
    
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

        # The chunks are given in the bam's order (shard order, then chunk order). The two
        # appends are independent copies (no recompression), so they run side by side.
        # `wait <pid>` returns that job's status, so a failed copy fails the task.
        #
        # Table: header, then each chunk's table as it is. Appended gzip members are one
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
            done < ~{write_lines(chunk_tsvs)}
        ) &
        tsv_pid=$!

        # Bam: the chunks cover disjoint reads, in the bam's order, so appending them
        # keeps the bam coordinate sorted; samtools cat copies the compressed blocks.
        samtools cat --no-PG -b ~{write_lines(chunk_bams)} -o ~{sample_id}.iso_cats.bam &
        bam_pid=$!

        wait $tsv_pid
        wait $bam_pid

        # Summed in chunk order, so categories stay in the order a single run first sees them.
        awk -F'\t' -v OFS='\t' '
            FNR == 1 { next }
            !($1 in count) { order[++n] = $1 }
            { count[$1] += $2 }
            END {
                print "Category", "Count"
                for (i = 1; i <= n; i++) print order[i], count[order[i]]
            }
        ' $(cat ~{write_lines(chunk_summaries)}) > ~{sample_id}.iso_cats.summary_counts.tsv

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
