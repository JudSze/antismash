from antismash.main import run_antismash, get_all_modules
from antismash.common.secmet import record

from antismash.config import build_config

from antismash.modules import nrps_pks
from antismash.modules.nrps_pks import specific_analysis

args = ["--minimal",
        "--enable-nrps-pks",
        "--genefinding-tool", "prodigal"]

options = build_config(args, isolated=True, modules=get_all_modules())

run_antismash("JF752342.1.gb", options)
gbk_record = record.Record.from_genbank("NC_003888.3.region021/DQ983361.1.region001.gbk")[0]
domain_res = nrps_pks.run_on_record(gbk_record, None, options)

chemical_structure = specific_analysis(gbk_record, domain_res, options)

