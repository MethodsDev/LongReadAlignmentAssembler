version 1.0

# Chunk-scatter LRAA: one genome-wide make-chunks task, one task per CHUNK, one
# global merge. This differs from LRAA.wdl, which scatters per CHROMOSOME and runs
# a complete `LRAA --chunk` inside each shard: here the unit of parallelism is the
# chunk itself, so a whole-genome run's makespan is one chunk's work rather than
# one chromosome's, and the leaves are small enough to run preemptible.
#
# STRANDLESS-ONLY on purpose: there is no chunk_by_strand input. Strandless is
# LRAA's default and the faster ordering, one chunk holds both orientations so
# prep writes half as many mini FASTAs and GTFs, and strand-first would require
# stage 1 -- a whole-BAM orientation split that would have to run in the prep
# task and be shipped to every leaf.
#
# LRAA.wdl stays as-is and remains the supported chromosome-scatter path.

workflow LRAA_chunk_scatter {
    input {
        String sample_id
        File referenceGenome
        File inputBAM
        File? inputBAMindex
        # Caller-supplied splice-graph evidence, already depth-normalized by the
        # caller. make_chunks partitions it per chunk EXACTLY as it partitions
        # inputBAM -- same region bounds, same mini contig, same offset, same
        # overhang rule -- so each leaf sees a splice-graph slice and a read slice
        # in one coordinate system. The slices land in the chunk directories and
        # ride to the leaves inside the per-chunk tars below.
        File? bam_for_sg
        File? bam_for_sg_index
        # ONE cut plan shared by SIBLING runs. Geometry only: which contig is cut
        # where. When supplied, make_chunks SKIPS cut selection and extracts on
        # this geometry, so every run handed the same plan produces identical
        # chunk bounds and therefore slices bam_for_sg identically.
        # num_total_reads and discovery are NEVER read from it -- both come from
        # THIS run, since the TPM denominator is per-caller.
        #
        # This is what makes bam_for_sg safe across callers. Cut POSITIONS are
        # otherwise per-caller: each run selects on its own --bam, gets its own
        # bounds, and slices the shared splice-graph evidence at those bounds, so
        # the runs are no longer comparable and their boundary overhang drops
        # differ. Selecting on bam_for_sg itself would be worse -- it is
        # depth-normalized, so a position that looks unspanned there can be
        # spanned by a raw read, which extraction then silently DROPS. ChunkedRun
        # refuses --bam_for_sg unless one of the shared-geometry inputs is given.
        File? chunk_plan
        File? annot_gtf
        Boolean quant_only = false
        Boolean HiFi = false
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        # Restrict the partition to these contigs, space-separated as
        # LRAA.wdl's main_chromosomes is. Empty means every reference the genome
        # fasta and the bam header agree on. Passed to ChunkedRun as --contigs,
        # which refuses a name absent from both rather than silently dropping it.
        String main_chromosomes = ""
        # Contigs to run in LRAA oversimplify mode, named as in the genome
        # fasta. ChunkedRun resolves it per chunk and rewrites the name to the
        # mini contig, so a caller writes the ORIGINAL name here.
        String? oversimplify
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        File? cell_list
        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
        Boolean stream_reads = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int? num_total_reads

        # make-chunks: one task running the genome-wide selection/extraction pool.
        # Its concurrency comes from --cpu_budget, so cpu IS that budget, and its
        # peak is therefore concurrency x per-unit rather than a constant.
        #
        # 8 GiB is PROVISIONAL, inferred rather than measured, and the distinction
        # matters here. On the PBMC whole-genome library (81.5M reads, 7.8 GB bam,
        # 305 chunks over 25 contigs) this task's timing.json reports PER-SUBPROCESS
        # peaks of 80 MiB for a cut selection and 90 MiB for an extraction, which at
        # 16 workers implies ~1.4 GiB. But that is the same interval-sampled,
        # one-subprocess-at-a-time metric that understated a LEAF's true container
        # peak by 4x (789 MB sampled against 3,003 MiB from the kernel), so the same
        # understatement is likely here and the implied 1.4 GiB is a floor, not a
        # peak. No task-level number exists for this task yet: the cgroup telemetry
        # that would give one (makeChunksResources) was added after that run.
        #
        # So 8 GiB would be a reduction resting on an inferred floor, not a
        # measured fit -- and this task is the wrong one to guess on: it is
        # non-preemptible, it holds the whole input, and every leaf and the merge
        # depend on it, so an under-request costs the entire run rather than one
        # retry. STAYS AT 32 GiB for the same reason the leaf default stays at 16:
        # both are reductions I can argue for and cannot yet verify, and this file
        # should not lower a production default on a number it documents as a lower
        # bound. 8 GiB is the candidate to apply once a whole-genome
        # makeChunksResources row exists -- the telemetry that produces it is now
        # wired, so this is one run away from being decidable.
        # Optional so a CALLER can forward its own knobs without restating these
        # numbers. LRAA.wdl passes Int? straight through; the values live here
        # once, resolved below.
        Int? makeChunksCpu
        Int? makeChunksMemoryGB

        # one task per chunk. See the note above the make_chunks call for why this
        # is a constant rather than derived from chunk width, and what the measured
        # peaks actually are.
        Int? chunkCpu
        Int? chunkMemoryGB
        Int chunkPreemptible = 3

        # stage 6, and now MEASURED at genome scale -- the one number here that
        # previously was not.
        #
        # Measured by running the merge directly over the 556 surviving units of a
        # whole-genome PBMC by_chunk run (278 chunks x 2 strands, 9,801,365 tracking
        # rows, 57,849 expr rows, coordinates translated): max RSS 371 MiB, wall
        # 2m53s.
        #
        # The number that makes this GENERALIZE rather than being one data point:
        # tracking_merge_peak_resident_rows came back at exactly 500,000, which is
        # the external sort's own run cap. The merge held AT the cap rather than
        # below it, so its peak is bounded by that cap and not by unit count or row
        # count -- a larger library spills more runs, it does not hold more rows.
        # That is the property the external sort was built for, and it is why a
        # constant is the right shape here where it was the wrong shape for a leaf.
        #
        # 8 GiB, not the ~4 the measurement alone would justify. This project's
        # rule for a memory reservation is a SAFE MINIMUM of 8 GiB wherever a
        # measurement suggests 4 or less: the asymmetry is not close. Over-asking
        # costs scheduling latency on a shared box and a little money on Terra;
        # under-asking a non-preemptible terminal gather discards every leaf's work
        # and the whole run's wall time. ~22x the measured 371 MiB, and still 4x
        # below the 32 GiB this replaced.
        #
        # Raise the external sort's run cap and this must move with it, since the
        # cap is what bounds the peak.
        Int? mergeCpu
        Int? mergeMemoryGB

        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    }

    # The defaults, in one place. Every number here is justified in the input
    # comments above. Prep and leaf memory are the PRE-EXISTING defaults, deliberately
    # unchanged: each is a reduction I can argue for from a lower-bound measurement
    # and cannot yet verify. The telemetry to settle them is wired.
    Int makeChunksCpu_use = select_first([makeChunksCpu, 16])
    Int makeChunksMemoryGB_use = select_first([makeChunksMemoryGB, 32])
    Int chunkCpu_use = select_first([chunkCpu, 2])
    Int chunkMemoryGB_use = select_first([chunkMemoryGB, 16])
    Int mergeCpu_use = select_first([mergeCpu, 2])
    Int mergeMemoryGB_use = select_first([mergeMemoryGB, 8])

    Boolean discovery = !quant_only
    String LRAA_output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    String LRAA_output_prefix = sample_id + "." + LRAA_output_suffix

    # NOTE on leaf memory, deliberately NOT derived from approx_MB_per_cut.
    #
    # MEASURED on the PBMC whole-genome library (81.5M reads, 305 chunks, 150
    # leaves sampled via cgroup memory.peak): median 1,037 MiB, p95 2,026, max
    # 3,003. So the 16 GiB below is ~5x the worst leaf of a real library at the
    # 10 Mb default, and on a ~300-leaf scatter that over-request is the
    # difference between a scatter that schedules and one that queues.
    #
    # A width-derived value (4 x (1 + 0.022 x span) GiB, mirroring
    # pylib/ChunkedRun.py's chunk_unit_peak_mib doubled for the two orientations)
    # lands at 5 GiB and fits the measurement -- but approx_MB_per_cut is a TARGET,
    # not a bound. The last segment of a contig absorbs the remainder, and a contig
    # whose every candidate position is blocked by annotation is not cut at all, so
    # an actual chunk can be far wider than the target: a whole chr1 unit measures
    # 5,594.7 MiB, 11.2 GiB for its pair, which a target-derived 5 GiB would
    # under-request precisely on the pathological chunk it has to cover.
    #
    # The correct fix is for make_chunks to emit a per-chunk memory value from each
    # chunk's ACTUAL span -- it already computes the spans, and chunk_unit_peak_mib
    # already exists -- as a file per chunk globbed in the same order as chunkTars,
    # so the scatter can zip them. Until that lands, the constant stays.

    call make_chunks {
        input:
            referenceGenome = referenceGenome,
            inputBAM = inputBAM,
            inputBAMindex = inputBAMindex,
            annot_gtf = annot_gtf,
            main_chromosomes = main_chromosomes,
            discovery = discovery,
            HiFi = HiFi,
            min_mapping_quality = min_mapping_quality,
            min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
            approx_MB_per_cut = approx_MB_per_cut,
            approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
            bam_for_sg = bam_for_sg,
            bam_for_sg_index = bam_for_sg_index,
            chunk_plan = chunk_plan,
            cell_list = cell_list,
            cell_barcode_tag = cell_barcode_tag,
            read_umi_tag = read_umi_tag,
            stream_reads = stream_reads,
            rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
            num_total_reads = num_total_reads,
            makeChunksCpu = makeChunksCpu_use,
            makeChunksMemoryGB = makeChunksMemoryGB_use,
            docker = docker
    }

    scatter (chunkTar in make_chunks.chunkTars) {
        call process_chunk {
            input:
                chunkTar = chunkTar,
                chunkPlan = make_chunks.chunkPlan,
                discovery = discovery,
                HiFi = HiFi,
                min_mapping_quality = min_mapping_quality,
                min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
                cell_list = cell_list,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                stream_reads = stream_reads,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                oversimplify = oversimplify,
                chunkCpu = chunkCpu_use,
                chunkMemoryGB = chunkMemoryGB_use,
                chunkPreemptible = chunkPreemptible,
                docker = docker
        }
    }

    call merge_chunks {
        input:
            unitsJsons = process_chunk.unitsJson,
            quantExprFiles = flatten(process_chunk.unitQuantExpr),
            quantTrackingFiles = flatten(process_chunk.unitQuantTracking),
            gtfFiles = flatten(process_chunk.unitGtf),
            readAssignmentSummaries = flatten(process_chunk.unitReadAssignmentSummary),
            discovery = discovery,
            outputFilePrefix = LRAA_output_prefix,
            mergeCpu = mergeCpu_use,
            mergeMemoryGB = mergeMemoryGB_use,
            docker = docker
    }

    output {
        File mergedQuantExpr = merge_chunks.mergedQuantExpr
        File mergedQuantTracking = merge_chunks.mergedQuantTracking
        File mergedTpmAudit = merge_chunks.mergedTpmAudit
        File mergedReadAssignmentSummary = merge_chunks.mergedReadAssignmentSummary
        File? mergedGTF = merge_chunks.mergedGTF
        File mergeResult = merge_chunks.mergeResult
        File chunkPlan = make_chunks.chunkPlan
        Array[File] chunkLogs = process_chunk.chunkLog
        Int numTotalReads = make_chunks.numTotalReads
        # What the run cost, as data rather than as something to reconstruct from
        # a log. One row per leaf plus the two singleton tasks, so a whole-genome
        # run answers "what should chunkMemoryGB be" directly. These are the
        # inputs to any future re-sizing; keep them delocalized.
        Array[File] chunkResources = process_chunk.chunkResources
        Array[File] unitResourceSummaries = flatten(process_chunk.unitResourceSummaries)
        File makeChunksResources = make_chunks.makeChunksResources
        File makeChunksTiming = make_chunks.makeChunksTiming
        File mergeResources = merge_chunks.mergeResources
    }
}


task make_chunks {
    input {
        File referenceGenome
        File inputBAM
        File? inputBAMindex
        File? bam_for_sg
        File? bam_for_sg_index
        File? chunk_plan
        File? annot_gtf
        Boolean discovery
        Boolean HiFi
        Int min_mapping_quality
        Int min_mapping_quality_for_final_quant
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        File? cell_list
        String cell_barcode_tag
        String read_umi_tag
        Boolean stream_reads
        Boolean rescue_unassigned_reads_via_transcriptome_alignment
        Int? num_total_reads
        String main_chromosomes
        Int makeChunksCpu
        Int makeChunksMemoryGB
        String docker
    }

    # 3x because the chunk directories hold a partition of the BAM plus its
    # indices alongside the localized input. The splice-graph bam is partitioned
    # the same way, so it counts on both sides of the multiplier.
    Float inputsGB = size(inputBAM, "GB") + size(referenceGenome, "GB") + (if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0) + (if defined(bam_for_sg) then size(bam_for_sg, "GB") else 0.0)
    Float diskRawGB = 3.0 * inputsGB + 50.0
    Int diskGB = if diskRawGB > 200.0 then ceil(diskRawGB) else 200
    String numTotalReadsStr = if defined(num_total_reads) then "~{select_first([num_total_reads])}" else ""

    command <<<
    set -euo pipefail

    # Symlink the inputs into the (writable) work dir: the input mounts can be
    # read-only, and both `samtools index` and faidx write beside their argument.
    mkdir -p inputs work
    ln -s ~{inputBAM} inputs/input.bam
    ~{if defined(inputBAMindex) then "ln -s " + select_first([inputBAMindex]) + " inputs/input.bam.bai" else ""}
    if [[ ! -e inputs/input.bam.bai && ! -e inputs/input.bam.csi ]]; then
        samtools index -@ ~{makeChunksCpu} inputs/input.bam
    fi
    ln -s ~{referenceGenome} inputs/genome.fa
    samtools faidx inputs/genome.fa
    ~{if defined(bam_for_sg) then "ln -s " + select_first([bam_for_sg]) + " inputs/sg.bam" else ""}
    ~{if defined(bam_for_sg_index) then "ln -s " + select_first([bam_for_sg_index]) + " inputs/sg.bam.bai" else ""}
    if [[ -e inputs/sg.bam && ! -e inputs/sg.bam.bai && ! -e inputs/sg.bam.csi ]]; then
        samtools index -@ ~{makeChunksCpu} inputs/sg.bam
    fi
    # Localized OUTSIDE work/: --stop_after_make_chunks writes this run's OWN
    # leaf plan at work/chunk_plan.json, which would clobber the shared plan
    # every sibling run reads.
    ~{if defined(chunk_plan) then "ln -s " + select_first([chunk_plan]) + " inputs/shared_cut_plan.json" else ""}

    # Same -F 0x904 policy as LRAA's own count_reads_from_bam, so a scattered
    # run's TPM denominator matches an unscattered one's. Computed here rather
    # than in its own task to avoid localizing the bam twice.
    N_OVERRIDE="~{numTotalReadsStr}"
    if [[ -n "${N_OVERRIDE}" ]]; then
        N="${N_OVERRIDE}"
    else
        N=$(samtools view -@ ~{makeChunksCpu} -c -F 0x904 inputs/input.bam)
    fi
    echo "num_total_reads=${N}" >&2
    echo "${N}" > work/num_total_reads.txt

    /usr/local/src/LRAA/pylib/ChunkedRun.py \
        --bam inputs/input.bam \
        ~{if defined(bam_for_sg) then "--bam_for_sg inputs/sg.bam" else ""} \
        ~{if defined(chunk_plan) then "--chunk_plan inputs/shared_cut_plan.json" else ""} \
        --genome_fa inputs/genome.fa \
        ~{if defined(annot_gtf) then "--gtf " + annot_gtf else ""} \
        ~{true="--discovery" false="" discovery} \
        --output_dir work \
        --cpu_budget ~{makeChunksCpu} \
        --num_total_reads "${N}" \
        --no_reuse_source_bam \
        --stop_after_make_chunks \
        ~{if main_chromosomes != "" then "--contigs " + sub(main_chromosomes, " +", ",") else ""} \
        ~{true="--HiFi" false="" HiFi} \
        --min_mapping_quality ~{min_mapping_quality} \
        --min_mapping_quality_for_final_quant ~{min_mapping_quality_for_final_quant} \
        ~{if defined(approx_MB_per_cut) then "--approx_MB_per_cut " + approx_MB_per_cut else ""} \
        ~{if defined(approx_MB_per_cut_wiggle_window) then "--approx_MB_per_cut_wiggle_window " + approx_MB_per_cut_wiggle_window else ""} \
        ~{true="--stream_reads" false="--no_stream_reads" stream_reads} \
        ~{true="" false="--no_stream_reads_rescue_unassigned" rescue_unassigned_reads_via_transcriptome_alignment} \
        ~{if defined(cell_list) then "--cell_list " + cell_list else ""} \
        ~{if (cell_barcode_tag != "CB") then "--cell_barcode_tag " + cell_barcode_tag else ""} \
        ~{if (read_umi_tag != "XM") then "--read_umi_tag " + read_umi_tag else ""}

    # One tar per chunk, so the scatter carries exactly one file per shard and every
    # chunk's identically-named members (chunk.bam, chunk.fa, ...) stay in their own
    # archive. Uncompressed: the mini bam is already compressed, and gzip here only
    # burns cpu. chunk_id is recovered downstream from the tar's basename.
    # The explicit member list is deliberate: it fails loudly on a chunk the
    # extractor did not fully write, rather than shipping a partial archive that a
    # leaf discovers hours later.
    mkdir -p staged
    python3 - <<'PY'
import json, os, subprocess
plan = json.load(open("work/chunk_plan.json"))
for entry in plan["chunks"]:
    cid = entry["chunk_id"]
    cdir = os.path.join("work", "chunks", cid)
    members = ["chunk.fa", "chunk.fa.fai", "chunk.bam", "chunk.bam.bai",
               "chunk.partition.json"]
    if entry["has_gtf"]:
        members.append("chunk.gtf")
    # The chunk's own manifest is the single source of truth for whether a
    # splice-graph slice exists as a LOCAL file: files.sg_bam is absent when no
    # --bam_for_sg was given, and under source-bam reuse it names the SOURCE bam
    # instead of a slice. --no_reuse_source_bam above means reuse cannot fire in
    # this shape, but keying on the manifest rather than on that flag keeps this
    # from being a second place that can disagree with the extractor.
    manifest = json.load(open(os.path.join(cdir, "chunk.partition.json")))
    if manifest["files"].get("sg_bam") and not manifest.get("sg_bam_reused_from_source"):
        members += ["chunk.sg.bam", "chunk.sg.bam.bai"]
    for m in members:
        p = os.path.join(cdir, m)
        if not os.path.exists(p):
            raise SystemExit("chunk {} is missing {}".format(cid, p))
    subprocess.run(["tar", "-cf", os.path.join("staged", cid + ".tar"),
                    "-C", cdir] + members, check=True)
print("tarred {} chunk(s)".format(len(plan["chunks"])))
PY

    # Same telemetry as the leaves, for the task whose peak is a FORMULA rather
    # than a constant: make_chunks runs its cut selections and extractions
    # concurrently under --cpu_budget, so its peak is roughly the per-unit peak
    # times the concurrency, and sizing it needs both numbers. Measured on the
    # PBMC whole-genome library (81.5M reads, 7.8 GB bam, 305 chunks over 25
    # contigs): per-unit peak 80 MiB for a cut selection and 90 MiB for an
    # extraction, so even at 16 workers this task's own high-water mark is far
    # below what it reserves.
    {
        printf 'task\tcontainer_peak_rss_bytes\tcpu\tmemory_gb_requested\tchunks\n'
        printf 'make_chunks\t%s\t%s\t%s\t%s\n' \
            "$(cat /sys/fs/cgroup/memory.peak 2>/dev/null || echo NA)" \
            "~{makeChunksCpu}" "~{makeChunksMemoryGB}" \
            "$(python3 -c 'import json;print(len(json.load(open("work/chunk_plan.json"))["chunks"]))' 2>/dev/null || echo NA)"
    } > make_chunks_resources.tsv
    >>>

    output {
        File chunkPlan = "work/chunk_plan.json"
        Array[File] chunkTars = glob("staged/*.tar")
        Int numTotalReads = read_int("work/num_total_reads.txt")
        # Its own timing.json carries the per-unit RSS and wall series this
        # summary's peak is the envelope of; both are wanted when re-sizing.
        File makeChunksResources = "make_chunks_resources.tsv"
        File makeChunksTiming = "work/timing.json"
    }

    runtime {
        docker: docker
        cpu: makeChunksCpu
        memory: "~{makeChunksMemoryGB} GiB"
        preemptible: 0
        bootDiskSizeGb: 30
        disks: "local-disk ~{diskGB} SSD"
    }
}


task process_chunk {
    input {
        File chunkTar
        File chunkPlan
        Boolean discovery
        Boolean HiFi
        Int min_mapping_quality
        Int min_mapping_quality_for_final_quant
        File? cell_list
        String cell_barcode_tag
        String read_umi_tag
        Boolean stream_reads
        Boolean rescue_unassigned_reads_via_transcriptome_alignment
        String? oversimplify
        Int chunkCpu
        Int chunkMemoryGB
        Int chunkPreemptible
        String docker
    }

    String chunkId = basename(chunkTar, ".tar")
    # 12x because stage 3b writes two strand BAMs, stage 4 two normalized BAMs,
    # and stage 5 its own scratch, all from the mini BAM.
    Float chunkDiskRawGB = 12.0 * size(chunkTar, "GB") + 20.0
    Int chunkDiskGB = if chunkDiskRawGB > 50.0 then ceil(chunkDiskRawGB) else 50

    command <<<
    set -euo pipefail

    mkdir -p work/chunks/~{chunkId} work/logs
    cp ~{chunkPlan} work/chunk_plan.json
    tar -xf ~{chunkTar} -C work/chunks/~{chunkId}

    /usr/local/src/LRAA/pylib/ChunkedRun.py \
        --output_dir work \
        --only_chunk ~{chunkId} \
        --cpu_budget ~{chunkCpu} \
        ~{true="--discovery" false="" discovery} \
        ~{true="--HiFi" false="" HiFi} \
        --min_mapping_quality ~{min_mapping_quality} \
        --min_mapping_quality_for_final_quant ~{min_mapping_quality_for_final_quant} \
        ~{true="--stream_reads" false="--no_stream_reads" stream_reads} \
        ~{true="" false="--no_stream_reads_rescue_unassigned" rescue_unassigned_reads_via_transcriptome_alignment} \
        ~{if defined(cell_list) then "--cell_list " + cell_list else ""} \
        ~{if (cell_barcode_tag != "CB") then "--cell_barcode_tag " + cell_barcode_tag else ""} \
        ~{if (read_umi_tag != "XM") then "--read_umi_tag " + read_umi_tag else ""}

    # Rename every unit artifact to <unit_id>.<suffix>. The merge addresses inputs as
    # quant_prefix + ".quant.expr" etc, so unit-scoped basenames let the merge task drop
    # every gathered file into ONE directory and set quant_prefix = <dir>/<unit_id>.
    # Without this, every chunk's files are named chunk_quant_plus.* and collide.
    mkdir -p staged
    python3 - <<'PY'
import json, os, shutil
doc = json.load(open(os.path.join("work", "chunks", "~{chunkId}", "units.json")))
for u in doc["units"]:
    uid, pfx = u["unit_id"], u["quant_prefix"]
    shutil.copy(pfx + ".quant.expr", os.path.join("staged", uid + ".quant.expr"))
    for suffix in (".quant.tracking.gz", ".quant.tracking"):
        if os.path.exists(pfx + suffix):
            shutil.copy(pfx + suffix,
                        os.path.join("staged", uid + ".quant.tracking.gz")
                        if suffix.endswith(".gz")
                        else os.path.join("staged", uid + ".quant.tracking"))
            break
    else:
        raise SystemExit("unit {} produced no quant.tracking".format(uid))
    if os.path.exists(pfx + ".gtf"):
        shutil.copy(pfx + ".gtf", os.path.join("staged", uid + ".gtf"))
    # UNCONDITIONAL, like quant.expr. Every work unit writes a read-assignment
    # summary now, including one with nothing to quantify -- run_quant_only's early
    # returns write an all-zero summary rather than returning before the writer. So
    # an absent file is a DEFECT, not a legitimate state, and stage 6 requires a
    # summary from every unit it merges. The `if os.path.exists` guard this replaces
    # was written against the old tolerance and is now dead code that would mask
    # exactly the missing input the merge refuses on purpose: the leaf would stage
    # nothing, the merge would report a unit with no summary, and the chunk that
    # actually failed to write one would be indistinguishable from a chunk that had
    # nothing to write.
    shutil.copy(pfx + ".read_assignment.summary.tsv",
                os.path.join("staged", uid + ".read_assignment.summary.tsv"))
shutil.copy(os.path.join("work", "chunks", "~{chunkId}", "units.json"),
            os.path.join("staged", "~{chunkId}.units.json"))
PY

    # What this task ACTUALLY cost, emitted as a delocalized output so it is
    # measurable on Terra rather than only in a local run tree. /sys/fs/cgroup/
    # memory.peak is the KERNEL's high-water mark for this container, so unlike
    # LRAA's own interval sampler it cannot miss a spike between samples --
    # measured on the PBMC whole-genome run, interval sampling reported 789 MB
    # where the cgroup peak for the same shape of task was 3,003 MiB, a 4x
    # understatement, because the sampler sees one unit at a time and the task
    # runs both orientations plus the driver. Sizing chunkMemoryGB off the
    # sampler would therefore have under-provisioned it.
    {
        printf 'chunk_id\tcontainer_peak_rss_bytes\tchunk_cpu\tchunk_memory_gb_requested\n'
        printf '%s\t%s\t%s\t%s\n' "~{chunkId}" \
            "$(cat /sys/fs/cgroup/memory.peak 2>/dev/null || echo NA)" \
            "~{chunkCpu}" "~{chunkMemoryGB}"
    } > "staged/~{chunkId}.chunk_resources.tsv"

    # The per-unit interval series LRAA writes for itself lives in the chunk work
    # directory, which is not delocalized. Copied out beside the cgroup peak: the
    # two answer different questions -- the peak sizes the request, the series
    # says which stage inside the unit spent it.
    for f in work/chunks/~{chunkId}/*.resources.summary.tsv; do
        [ -e "$f" ] || continue
        cp "$f" "staged/~{chunkId}.$(basename "$f")"
    done
    >>>

    output {
        File unitsJson = "staged/~{chunkId}.units.json"
        Array[File] unitQuantExpr = glob("staged/*.quant.expr")
        Array[File] unitQuantTracking = glob("staged/*.quant.tracking*")
        Array[File] unitGtf = glob("staged/*.gtf")
        Array[File] unitReadAssignmentSummary = glob("staged/*.read_assignment.summary.tsv")
        File chunkLog = "work/logs/chunk_~{chunkId}.log"
        # Resource telemetry, one row per chunk plus LRAA's own per-unit series.
        # Gathered by the caller so a whole-genome run yields a table of what its
        # ~300 leaves cost, which is what chunkMemoryGB should be set from.
        File chunkResources = "staged/~{chunkId}.chunk_resources.tsv"
        Array[File] unitResourceSummaries = glob("staged/*.resources.summary.tsv")
    }

    runtime {
        docker: docker
        cpu: chunkCpu
        memory: "~{chunkMemoryGB} GiB"
        preemptible: chunkPreemptible
        disks: "local-disk ~{chunkDiskGB} HDD"
    }
}


task merge_chunks {
    input {
        Array[File] unitsJsons
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        Array[File] gtfFiles
        Array[File] readAssignmentSummaries
        Boolean discovery
        String outputFilePrefix
        Int mergeCpu
        Int mergeMemoryGB
        String docker
    }

    Float mergeInputsGB = size(quantExprFiles, "GB") + size(quantTrackingFiles, "GB") + size(gtfFiles, "GB")
    Float mergeDiskRawGB = 3.0 * mergeInputsGB + 20.0
    Int mergeDiskGB = if mergeDiskRawGB > 100.0 then ceil(mergeDiskRawGB) else 100

    command <<<
    set -euo pipefail

    mkdir -p staged work
    for f in ~{sep=' ' quantExprFiles} ~{sep=' ' quantTrackingFiles} ~{sep=' ' gtfFiles} ~{sep=' ' readAssignmentSummaries}; do
        cp "$f" "staged/$(basename "$f")"
    done

    # Manifest order must reproduce ChunkedRun.ordered_units exactly:
    # strand + before -, then each unit's own `order`. Order is load-bearing for
    # the merge. The duplicate-unit_id check exists because a retried preemptible
    # shard whose outputs were gathered twice would otherwise double-count reads
    # in the merged tables with no error -- the merge validates offsets but not
    # uniqueness.
    python3 - <<'PY'
import json, os
paths = """~{sep="\n" unitsJsons}""".split()
units = []
for p in paths:
    for u in json.load(open(p))["units"]:
        units.append(u)
rank = {"+": 0, "-": 1}
units.sort(key=lambda u: (rank[u["strand"]], u["order"]))
seen = set()
out = []
for u in units:
    uid = u["unit_id"]
    if uid in seen:
        raise SystemExit("duplicate unit_id {} across chunk shards".format(uid))
    seen.add(uid)
    expr = os.path.abspath(os.path.join("staged", uid + ".quant.expr"))
    if not os.path.exists(expr):
        raise SystemExit("unit {} has no staged {}".format(uid, expr))
    out.append({
        "unit_id": uid,
        "quant_prefix": os.path.abspath(os.path.join("staged", uid)),
        "offset": int(u["offset"]),
    })
json.dump({"units": out}, open("manifest.json", "wt"), indent=2)
print("manifest: {} unit(s)".format(len(out)))
PY

    merge_chunk_outputs.py \
        --manifest manifest.json \
        --outdir work \
        ~{true="--discovery" false="" discovery} \
        --result "~{outputFilePrefix}.merge_result.json"

    mv work/merged/chunked.quant.expr "~{outputFilePrefix}.quant.expr"
    mv work/merged/chunked.quant.tracking.gz "~{outputFilePrefix}.quant.tracking.gz"
    mv work/merged/chunked.quant.expr.tpm_chunk_local.tsv "~{outputFilePrefix}.quant.expr.tpm_chunk_local.tsv"
    mv work/merged/chunked.read_assignment.summary.tsv "~{outputFilePrefix}.read_assignment.summary.tsv"
    if [[ "~{discovery}" == "true" ]]; then
        mv work/merged/chunked.gtf "~{outputFilePrefix}.gtf"
    fi

    # The merge is where a whole-genome run's peak lives: it pools every unit's
    # tracking rows through an external sort, and the number that decides its
    # ceiling is tracking_merge_peak_resident_rows in the merge result. Emitting
    # the container's own high-water mark beside the row count gives both halves
    # of the sizing -- bytes per resident row is what makes mergeMemoryGB
    # predictable for a library of any size, rather than a constant chosen once.
    {
        printf 'task\tcontainer_peak_rss_bytes\tcpu\tmemory_gb_requested\tunits\tpeak_resident_rows\n'
        printf 'merge_chunks\t%s\t%s\t%s\t%s\t%s\n' \
            "$(cat /sys/fs/cgroup/memory.peak 2>/dev/null || echo NA)" \
            "~{mergeCpu}" "~{mergeMemoryGB}" \
            "$(python3 -c 'import json;print(len(json.load(open("manifest.json"))["units"]))' 2>/dev/null || echo NA)" \
            "$(python3 -c 'import json,sys;d=json.load(open("~{outputFilePrefix}.merge_result.json"));print(d.get("tracking_merge_peak_resident_rows","NA"))' 2>/dev/null || echo NA)"
    } > merge_chunks_resources.tsv
    >>>

    output {
        File mergedQuantExpr = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTracking = "~{outputFilePrefix}.quant.tracking.gz"
        File mergedTpmAudit = "~{outputFilePrefix}.quant.expr.tpm_chunk_local.tsv"
        File mergedReadAssignmentSummary = "~{outputFilePrefix}.read_assignment.summary.tsv"
        File? mergedGTF = "~{outputFilePrefix}.gtf"
        File mergeResult = "~{outputFilePrefix}.merge_result.json"
        File mergeResources = "merge_chunks_resources.tsv"
    }

    runtime {
        docker: docker
        cpu: mergeCpu
        memory: "~{mergeMemoryGB} GiB"
        preemptible: 0
        disks: "local-disk ~{mergeDiskGB} SSD"
    }
}
