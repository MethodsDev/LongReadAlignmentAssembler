#!/usr/bin/env python3

"""Every Terra task names a machine that holds its request, next to cpu and memory.

Rules from the Batch probes of 2026-10-30 and from miniwdl:
- predefinedMachineType may sit beside cpu/memory, but Batch rejects a request larger
  than the machine, and rejects 8 GiB on c3d-highcpu-4 although the shape lists 8 GiB;
- miniwdl ignores predefinedMachineType and reserves 0 cores and 0 memory for a task with
  no cpu/memory, so those lines must stay;
- an HDD task cannot go on C3D, which Compute Engine does not list with pd-standard (not
  probed; see docs/v0.47.1_terra_runtime_plan.md).
The tier text is repeated per task because WDL 1.0 has no functions, so it is checked here.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
WDL_FILES = [
    "WDL/LRAA.wdl",
    "WDL/LRAA-cell_cluster_guided.wdl",
    "WDL/LRAA_quant_by_cluster.wdl",
    "WDL/subwdls/Incorporate_gene_symbols.wdl",
    "WDL/subwdls/LRAA-build_sparse_matrices_from_tracking.wdl",
    "WDL/subwdls/LRAA-filter_good_cells.wdl",
    "WDL/subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl",
    "WDL/subwdls/LRAA_chunk_scatter.wdl",
    "WDL/subwdls/LRAA_runner.wdl",
    "WDL/subwdls/Normalize_bam.wdl",
    "WDL/subwdls/Partition_data_by_chromosome.wdl",
    "WDL/subwdls/partition_bam_by_cell_cluster.wdl",
]

# guestCpus and GiB from the Compute Engine machineTypes list, us-central1-a, 2026-10-31
C3D_CORES = (4, 8, 16, 30, 60, 90, 180)
C3D_RAM = {"standard": lambda c: c * 4, "highmem": lambda c: c * 8}
C3D_HIGHCPU_RAM = {4: 8, 8: 16, 16: 32, 30: 59, 60: 118, 90: 177, 180: 354}
N2D_CORES = (2, 4, 8, 16, 32, 48, 64, 80, 96)
N2D_RAM = {"highcpu": 1, "standard": 4, "highmem": 8}


def _capacity(machine):
    family, kind, cores = machine.split("-")
    cores = int(cores)
    if family == "c3d":
        assert cores in C3D_CORES, machine
        return cores, C3D_HIGHCPU_RAM[cores] if kind == "highcpu" else C3D_RAM[kind](cores)
    assert family == "n2d" and cores in N2D_CORES, machine
    return cores, cores * N2D_RAM[kind]


def _task_bodies():
    for rel in WDL_FILES:
        src = (REPO / rel).read_text()
        heads = list(re.finditer(r"^(task|workflow)\s+(\w+)", src, re.M))
        for i, head in enumerate(heads):
            if head.group(1) == "task":
                end = heads[i + 1].start() if i + 1 < len(heads) else len(src)
                yield rel, head.group(2), src[head.start() : end]


TASKS = list(_task_bodies())
IDS = [name for _, name, _ in TASKS]


def _runtime(body):
    return re.search(r"^\s*runtime\s*\{(.*?)^\s*\}", body, re.M | re.S).group(1)


def _shape_block(body):
    start = re.search(r"^    # Machine shape \(v0\.47\.1", body, re.M).start()
    return body[start : re.search(r"^\s*command\s*<<<", body, re.M).start()]


@pytest.mark.parametrize("rel,name,body", TASKS, ids=IDS)
def test_machine_sits_beside_the_reservation(rel, name, body):
    runtime = _runtime(body)
    assert "predefinedMachineType: machine_type" in runtime
    assert re.search(r"^\s*memory:", runtime, re.M)
    # run_filter_good_cells has never declared cpu; it is not a reservation this release adds
    if name != "run_filter_good_cells":
        assert re.search(r"^\s*cpu:", runtime, re.M)


@pytest.mark.parametrize("rel,name,body", TASKS, ids=IDS)
def test_an_hdd_task_is_not_on_c3d(rel, name, body):
    if '" HDD"' in _runtime(body) or re.search(r"\sHDD\"", _runtime(body)):
        assert '"c3d-' not in _shape_block(body)


def test_the_shape_text_is_identical_across_tasks():
    variants = set()
    for _, _, body in TASKS:
        block = _shape_block(body)
        block = re.sub(r"Int shape_cpu = .*", "Int shape_cpu = X", block)
        block = re.sub(r"Int shape_memGB = .*", "Int shape_memGB = X", block)
        variants.add(block)
    assert len(variants) == 2  # one C3D text, one N2D text


def _eval_machine(task, cpu, mem):
    WDL = pytest.importorskip("WDL")

    env = WDL.Env.Bindings()
    env = env.bind("shape_cpu", WDL.Value.Int(cpu)).bind("shape_memGB", WDL.Value.Int(mem))
    stdlib = WDL.StdLib.Base("1.0")
    for decl in task.postinputs:
        if decl.name.startswith("shape_") and decl.name not in ("shape_cpu", "shape_memGB"):
            env = env.bind(decl.name, decl.expr.eval(env, stdlib))
    machine = [d for d in task.postinputs if d.name == "machine_type"][0]
    return machine.expr.eval(env, stdlib).value


@pytest.mark.parametrize(
    "rel,task_name", [("WDL/LRAA.wdl", "count_bam"), ("WDL/LRAA.wdl", "derive_contigs")]
)
def test_the_machine_always_holds_the_request(rel, task_name):
    """count_bam carries the C3D text and derive_contigs the N2D text."""

    WDL = pytest.importorskip("WDL")
    task = [t for t in WDL.load(str(REPO / rel)).tasks if t.name == task_name][0]
    top_cpu, top_mem = (180, 1440) if task_name == "count_bam" else (96, 768)
    for cpu in (1, 2, 3, 4, 5, 8, 9, 16, 17, 30, 32, 48, 60, 64, 80, 90, 96):
        for mem in (1, 2, 4, 7, 8, 9, 15, 16, 17, 24, 32, 33, 59, 60, 64, 65, 128, 200, 240, 256, 384, 480, 640, 720):
            if cpu > top_cpu or mem > top_mem:
                continue
            machine = _eval_machine(task, cpu, mem)
            cores, ram = _capacity(machine)
            assert cores >= cpu and ram >= mem, (cpu, mem, machine)
            if machine == "c3d-highcpu-4":
                assert mem <= 7, (cpu, mem)


def test_the_batch_rejection_cases():
    WDL = pytest.importorskip("WDL")
    tasks = {t.name: t for t in WDL.load(str(REPO / "WDL/LRAA.wdl")).tasks}
    assert _eval_machine(tasks["count_bam"], 4, 8) == "c3d-standard-4"
    assert _eval_machine(tasks["count_bam"], 4, 7) == "c3d-highcpu-4"
    assert _eval_machine(tasks["count_bam"], 16, 32) == "c3d-highcpu-16"
    assert _eval_machine(tasks["derive_contigs"], 1, 2) == "n2d-highcpu-2"
    assert _eval_machine(tasks["derive_contigs"], 1, 16) == "n2d-highmem-2"


def test_runner_is_spot_and_the_terminal_tasks_are_not():
    bodies = {name: body for _, name, body in TASKS}
    assert "Int preemptible_tries = 2" in bodies["LRAA_runner_task"]
    assert "preemptible: preemptible_tries" in _runtime(bodies["LRAA_runner_task"])
    for name in ("make_chunks", "merge_chunks", "emit_shared_chunk_plan"):
        assert "preemptible: 0" in _runtime(bodies[name]), name
