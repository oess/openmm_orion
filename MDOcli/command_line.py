import click

from openeye import oechem

from datarecord import read_mol_record, OEField

from MDOrion.Standards import Fields

from datarecord import (Types,
                        OEWriteRecord,
                        OERecord)

from datarecord.datarecord import RecordVecData, RecordData

from orionclient.types import File, Shard, ShardCollection

import os

from orionclient.session import (OrionSession,
                                 get_profile_config,
                                 get_session,
                                 APISession)

from MDOrion.Standards.mdrecord import MDDataRecord

from tempfile import TemporaryDirectory

import pickle

import parmed

@click.group(
    context_settings={
        "help_option_names": ("-h", "--help")
    }
)
@click.pass_context
def main(ctx):
    ctx.obj = dict()


@main.group()
@click.argument('filename', type=click.Path(exists=True))
@click.option("--id", help="Record ID number", default="all")
@click.option("--profile", help="OCLI profile name", default="default")
@click.pass_context
def dataset(ctx, filename, id, profile=None, max_retries=2):
    """Records Extraction"""

    ctx.obj['filename'] = filename

    ifs = oechem.oeifstream(filename)

    records = []

    while True:
        record = read_mol_record(ifs)
        if record is None:
            break
        records.append(record)
    ifs.close()

    if id == 'all':
        ctx.obj['records'] = records
    else:
        if int(id) < len(records):
            ctx.obj['records'] = [records[int(id)]]
        else:
            raise ValueError("Wrong record number selection: {} > max = {}".format(int(id), len(records)))
    # TODO
    if profile == "default" and os.environ.get("ORION_PROFILE") is not None:
        profile = os.environ["ORION_PROFILE"]

    profile_config = get_profile_config(profile=profile)

    ctx.obj['profile'] = profile
    ctx.obj['session'] = OrionSession(
            config=profile_config,
            requests_session=get_session({404: max_retries}))

    ctx.obj['credentials'] = profile_config

    if ctx.obj["credentials"] is not None:
        ctx.obj['session'].config = ctx.obj["credentials"]
    else:
        click.secho("Unable to find credentials", fg='red', err=True)

#
# @dataset.command("trajectory")
# @click.option("--format", help="Trajectory format", default='tar.gz')
# @click.option("--stgn", help="MD Stage number", default="all")
# @click.option("--fixname", help="Edit the trajectory file name", default=None)
# @click.pass_context
# def trajectory_extraction(ctx, format, stgn, fixname):
#     stagen = stgn
#
#     new_record_list = []
#
#     for record in ctx.obj['records']:
#
#         if not record.has_value(Fields.md_stages):
#             print("No MD stages have been found in the selected record")
#             continue
#
#         stages = record.get_value(Fields.md_stages)
#         nstages = len(stages)
#         title = record.get_value(Fields.title)
#         id = record.get_value(Fields.id)
#         fn = title + "_" + str(id)
#
#         if stagen == 'last':
#             stages_work_on = [stages[-1]]
#         elif stagen == 'all':
#             stages_work_on = stages
#         else:
#             if int(stagen) > nstages:
#                 print("Wrong stage number selection: {} > max = {}".format(stagen, nstages))
#                 return
#             else:
#                 stages_work_on = [stages[int(stagen)]]
#
#         for idx in range(0, len(stages_work_on)):
#
#             stage = stages_work_on[idx]
#
#             if stage.has_value(Fields.trajectory):
#                 fnt = stage.get_value(Fields.trajectory)
#                 trj_field = stage.get_field(Fields.trajectory.get_name())
#                 trj_meta = trj_field.get_meta()
#
#             elif stage.has_value(Fields.orion_local_trj_field):
#                 trj_id = stage.get_value(Fields.orion_local_trj_field)
#                 trj_field = stage.get_field(Fields.orion_local_trj_field.get_name())
#                 trj_meta = trj_field.get_meta()
#
#                 suffix = ''
#                 if stage.has_value(Fields.log_data):
#                     log = stage.get_value(Fields.log_data)
#                     log_split = log.split()
#                     for i in range(0, len(log_split)):
#                         if log_split[i] == 'suffix':
#                             suffix = log_split[i+2]
#                             break
#                 fnt = fn + '-' + suffix + '.' + format
#
#                 resource = ctx.obj['session'].get_resource(File, trj_id)
#
#                 resource.download_to_file(fnt)
#
#             else:
#                 print("No MD trajectory found in the selected stage record {}".format(stage.get_value(Fields.stage_name)))
#                 continue
#
#             if fixname is not None:
#                 trj_field = OEField(Fields.trajectory.get_name(),
#                                     Fields.trajectory.get_type(),
#                                     meta=trj_meta)
#
#                 stage.set_value(trj_field, fnt)
#
#                 stages_work_on[idx] = stage
#
#             record.set_value(Fields.md_stages, stages_work_on)
#             new_record_list.append(record)
#
#     if fixname is not None:
#
#         ofs = oechem.oeofstream(fixname)
#
#         for record in new_record_list:
#             OEWriteRecord(ofs, record, fmt='binary')
#


@dataset.command("makelocal")
@click.option("--name", help="Edit the trajectory file name", default="local.oedb")
@click.pass_context
def data_trajectory_extraction(ctx, name):

    new_records = []

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        new_record = OERecord(record)

        if not record.has_field(Fields.collection):
            raise ValueError("No Collection field has been found in the record")

        session = APISession

        collection_id = record.get_value(Fields.collection)

        collection = session.get_resource(ShardCollection, collection_id)

        stages = mdrecord.get_stages

        system_title = mdrecord.get_title
        sys_id = mdrecord.get_id

        new_stages = []

        for stage in stages:

            stg_type = stage.get_value(Fields.stage_type)
            new_stage = OERecord(stage)

            with TemporaryDirectory() as output_directory:
                data_fn = os.path.basename(output_directory) + '_' + system_title + '_' + str(sys_id) + '-' + stg_type + '.tar.gz'
                shard_id = stage.get_value(OEField("MDData_OPLMD", Types.Int))
                shard = session.get_resource(Shard(collection=collection), shard_id)
                shard.download_to_file(data_fn)
                new_stage.set_value(Fields.mddata, data_fn)

                if stage.has_field(OEField("Trajectory_OPLMD", Types.Int)):
                    trj_id = stage.get_value(OEField("Trajectory_OPLMD", Types.Int))
                    trj_fn = os.path.basename(output_directory) + '_' + system_title + '_' + str(sys_id) + '-' + stg_type + '_traj' + '.tar.gz'
                    resource = session.get_resource(File, trj_id)
                    resource.download_to_file(trj_fn)
                    new_stage.set_value(Fields.trajectory, trj_fn)

            new_stages.append(new_stage)

        new_record.set_value(Fields.md_stages, new_stages)

        if record.has_field(OEField('Structure_Parmed_OPLMD', Types.Int)):
            pmd_id = record.get_value(OEField('Structure_Parmed_OPLMD', Types.Int))
            shard = session.get_resource(Shard(collection=collection), pmd_id)

            with TemporaryDirectory() as output_directory:
                parmed_fn = os.path.join(output_directory, "parmed.pickle")

                shard.download_to_file(parmed_fn)

                with open(parmed_fn, 'rb') as f:
                    parm_dic = pickle.load(f)

                pmd_structure = parmed.structure.Structure()
                pmd_structure.__setstate__(parm_dic)

            new_record.set_value(Fields.pmd_structure, pmd_structure)

        if record.has_field(OEField('OETraj', Types.Record)):

            oetrajrec = record.get_value(OEField('OETraj', Types.Record))

            prot_conf_id = oetrajrec.get_value(OEField("ProtTraj_OPLMD", Types.Int))

            shard = session.get_resource(Shard(collection=collection),  prot_conf_id)

            with TemporaryDirectory() as output_directory:
                protein_fn = os.path.join(output_directory, "prot_traj_confs.oeb")

                shard.download_to_file(protein_fn)

                protein_conf = oechem.OEMol()

                with oechem.oemolistream(protein_fn) as ifs:
                    oechem.OEReadMolecule(ifs, protein_conf)

            oetrajrec.set_value(Fields.protein_traj_confs, protein_conf)

            new_record.set_value(OEField('OETraj', Types.Record), oetrajrec)

        new_record.delete_field(Fields.collection)

        new_records.append(new_record)

    ofs = oechem.oeofstream(name)

    for rec in new_records:
        OEWriteRecord(ofs, rec, fmt='binary')


@dataset.command("logs")
@click.option("--stgn", help="MD Stage name", default="last")
@click.pass_context
def logs_extraction(ctx, stgn):

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        info = mdrecord.get_stage_info(stg_name=stgn)

        print(info)


@dataset.command("protein")
@click.pass_context
def protein_extraction(ctx):

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        if not record.has_value(Fields.protein):
            print("No protein have been found in the selected record")
            return
        else:
            title = mdrecord.get_title
            fn = title.split('_')[0] + ".oeb"
            with oechem.oemolostream(fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, record.get_value(Fields.protein))
        print("Protein file generated: {}".format(fn))


@dataset.command("ligand")
@click.pass_context
def ligand_extraction(ctx):

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        if not record.has_value(Fields.ligand):
            print("No ligand have been found in the selected record")
            return
        else:
            title = mdrecord.get_title.split("_")[1:]
            title = "_".join(title)
            id = mdrecord.get_id
            fn = title + "_" + str(id)+".oeb"
            with oechem.oemolostream(fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, record.get_value(Fields.ligand))
        print("Ligand file generated: {}".format(fn))


@dataset.command("info")
@click.pass_context
def info_extraction(ctx):

    def GetHumanReadable(size, precision=2):
        suffixes = ['B', 'KB', 'MB', 'GB', 'TB']
        suffixIndex = 0

        while size > 1024 and suffixIndex < 4:
            suffixIndex += 1  # increment the index of the suffix
            size = size / 1024.0  # apply the division
        return "%.*f %s" % (precision, size, suffixes[suffixIndex])

    def recursive_record(record, level=0):

        for field in record.get_fields():

            field_type = field.get_type()

            blank = "       "
            print("{} |".format(blank * (level + 1)))
            print("{} |".format(blank * (level + 1)))
            dis = "______"

            if not field_type == RecordData and not field_type == RecordVecData:

                if (field.get_type() is Types.String or
                        field.get_type() is Types.Int or
                        field.get_type() is Types.Float):
                    print("{} {} name = {}\n        "
                          "{}type = {}\n        "
                          "{}value = {}\n        "
                          "{}size = {}".format(blank * (level + 1),
                                               dis,
                                               field.get_name(),
                                               blank * (level + 1),
                                               field.get_type(),
                                               blank * (level + 1),
                                               str(record.get_value(field))[0:30],
                                               blank * (level + 1),
                                               GetHumanReadable(record.get_value_size(field))
                                               ))
                else:
                    print("{} {} name = {}\n        "
                          "{}type = {}\n        "
                          "{}size = {}".format(blank * (level + 1),
                                               dis,
                                               field.get_name(),
                                               blank * (level + 1),
                                               field.get_type(),
                                               blank * (level + 1),
                                               GetHumanReadable(record.get_value_size(field))
                                               ))

            elif field_type == RecordData:
                print("{} {} RECORD: {}".format(blank * (level + 1), dis, field.get_name()))
                recursive_record(record.get_value(field), level + 1)

            elif field_type == RecordVecData:
                vec = record.get_value(field)
                print("{} {} RECORD VECTOR: {} containing {} records".format(blank * (level + 1),
                                                                             dis,
                                                                             field.get_name(),
                                                                             len(vec)))
                print("{} |".format(blank * (level + 2)))
                print("{} |".format(blank * (level + 2)))

                for idx in range(0, len(vec)):
                    print("{} {} RECORD # {}".format(blank * (level + 2),
                                                     dis,
                                                     idx))

                    recursive_record(vec[idx], level + 2)
                    if idx != len(vec) - 1:
                        print("{} |".format(blank * (level + 2)))
                        print("{} |".format(blank * (level + 2)))
            else:
                raise ValueError("Field type error: {}".format(field_type))

    for idx in range(0, len(ctx.obj['records'])):
        print(30 * "*" + " RECORD {}/{} ".format(idx + 1, len(ctx.obj['records'])) + 30 * "*")
        recursive_record(ctx.obj['records'][idx], 0)
        print("\n" + 30 * "*" + " END RECORD ".format(idx + 1, len(ctx.obj['records'])) + 30 * "*" + "\n") 