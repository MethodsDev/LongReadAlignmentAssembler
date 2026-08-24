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

workflow LRAA_chunk_scatter_wf {
    input {
        String sample_id
        File referenceGenome
        File inputBAM
        File? inputBAMindex
        File? annot_gtf
        Boolean quant_only = false
        Boolean HiFi = false
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        File? cell_list
        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
        Boolean stream_reads = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int? num_total_reads

        # make-chunks: one task running the genome-wide selection/extraction pool.
        # Its concurrency comes from --cpu_budget, so cpu IS that budget.
        Int makeChunksCpu = 16
        Int makeChunksMemoryGB = 32

        # one task per chunk. 16 GiB is the requested box, not a measured peak:
        # ChunkedRun's own scheduler ESTIMATES a chunk at
        # CHUNK_UNIT_FIXED_MIB + 22 MiB per genomic Mb (pylib/ChunkedRun.py),
        # which is a deliberate upper envelope, so treat 16 as provisioned headroom
        # and raise it if a leaf is OOM-killed.
        Int chunkCpu = 2
        Int chunkMemoryGB = 16
        Int chunkPreemptible = 3

        # stage 6. Sized above the leaves because the merge is where a whole-genome
        # run's peak RSS lives -- an external merge sort over every unit's tracking rows.
        Int mergeCpu = 2
        Int mergeMemoryGB = 32

        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    }

    Boolean discovery = !quant_only
    String LRAA_output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    String LRAA_output_prefix = sample_id + "." + LRAA_output_suffix

    call make_chunks {
        input:
            referenceGenome = referenceGenome,
            inputBAM = inputBAM,
            inputBAMindex = inputBAMindex,
            annot_gtf = annot_gtf,
            discovery = discovery,
            HiFi = HiFi,
            min_mapping_quality = min_mapping_quality,
            min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
            approx_MB_per_cut = approx_MB_per_cut,
            approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
            cell_list = cell_list,
            cell_barcode_tag = cell_barcode_tag,
            read_umi_tag = read_umi_tag,
            stream_reads = stream_reads,
            rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
            num_total_reads = num_total_reads,
            makeChunksCpu = makeChunksCpu,
            makeChunksMemoryGB = makeChunksMemoryGB,
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
                chunkCpu = chunkCpu,
                chunkMemoryGB = chunkMemoryGB,
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
            discovery = discovery,
            outputFilePrefix = LRAA_output_prefix,
            mergeCpu = mergeCpu,
            mergeMemoryGB = mergeMemoryGB,
            docker = docker
    }

    output {
        File mergedQuantExpr = merge_chunks.mergedQuantExpr
        File mergedQuantTracking = merge_chunks.mergedQuantTracking
        File mergedTpmAudit = merge_chunks.mergedTpmAudit
        File? mergedGTF = merge_chunks.mergedGTF
        File mergeResult = merge_chunks.mergeResult
        File chunkPlan = make_chunks.chunkPlan
        Array[File] chunkLogs = process_chunk.chunkLog
        Int numTotalReads = make_chunks.numTotalReads
    }
}


task make_chunks {
    input {
        File referenceGenome
        File inputBAM
        File? inputBAMindex
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
        Int makeChunksCpu
        Int makeChunksMemoryGB
        String docker
    }

    # 3x because the chunk directories hold a partition of the BAM plus its
    # indices alongside the localized input.
    Float inputsGB = size(inputBAM, "GB") + size(referenceGenome, "GB") + (if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0)
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
        --genome_fa inputs/genome.fa \
        ~{if defined(annot_gtf) then "--gtf " + annot_gtf else ""} \
        ~{true="--discovery" false="" discovery} \
        --output_dir work \
        --cpu_budget ~{makeChunksCpu} \
        --num_total_reads "${N}" \
        --no_reuse_source_bam \
        --stop_after_make_chunks \
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
    for m in members:
        p = os.path.join(cdir, m)
        if not os.path.exists(p):
            raise SystemExit("chunk {} is missing {}".format(cid, p))
    subprocess.run(["tar", "-cf", os.path.join("staged", cid + ".tar"),
                    "-C", cdir] + members, check=True)
print("tarred {} chunk(s)".format(len(plan["chunks"])))
PY
    >>>

    output {
        File chunkPlan = "work/chunk_plan.json"
        Array[File] chunkTars = glob("staged/*.tar")
        Int numTotalReads = read_int("work/num_total_reads.txt")
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
shutil.copy(os.path.join("work", "chunks", "~{chunkId}", "units.json"),
            os.path.join("staged", "~{chunkId}.units.json"))
PY
    >>>

    output {
        File unitsJson = "staged/~{chunkId}.units.json"
        Array[File] unitQuantExpr = glob("staged/*.quant.expr")
        Array[File] unitQuantTracking = glob("staged/*.quant.tracking*")
        Array[File] unitGtf = glob("staged/*.gtf")
        File chunkLog = "work/logs/chunk_~{chunkId}.log"
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
    for f in ~{sep=' ' quantExprFiles} ~{sep=' ' quantTrackingFiles} ~{sep=' ' gtfFiles}; do
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
    if [[ "~{discovery}" == "true" ]]; then
        mv work/merged/chunked.gtf "~{outputFilePrefix}.gtf"
    fi
    >>>

    output {
        File mergedQuantExpr = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTracking = "~{outputFilePrefix}.quant.tracking.gz"
        File mergedTpmAudit = "~{outputFilePrefix}.quant.expr.tpm_chunk_local.tsv"
        File? mergedGTF = "~{outputFilePrefix}.gtf"
        File mergeResult = "~{outputFilePrefix}.merge_result.json"
    }

    runtime {
        docker: docker
        cpu: mergeCpu
        memory: "~{mergeMemoryGB} GiB"
        preemptible: 0
        disks: "local-disk ~{mergeDiskGB} SSD"
    }
}
