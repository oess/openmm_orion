from datarecord import Types, OEField, OEFieldMeta, Meta

from orionclient.session import in_orion

# Orion Hidden meta data options
_metaHidden = OEFieldMeta(options=[Meta.Display.Hidden])
_metaIDHidden = OEFieldMeta(options=[Meta.Source.ID, Meta.Display.Hidden])


class Fields:
    # Current number of MD steps
    current_iteration_field = OEField("Current_Iterations_OMD", Types.Int)

    # Total number of MD steps
    md_nsteps_field = OEField("MD_nsteps_OMD", Types.Int)

    # Current number of cycles
    cycle_id = OEField("Cycle_ID_OMD", Types.Int)

    # Tpr binary file
    tpr_field = OEField("TPR_bytes_OMD", Types.Blob, meta=_metaHidden)

    # Prefix name field
    prefix_name_field = OEField("Prefix_OPLMD", Types.String)

    if in_orion():
        trajectory = OEField("GMXTrajectory_OMD", Types.Int, meta=_metaHidden)
        gmx_restart = OEField("GMXRestart_OMD", Types.Int, meta=_metaHidden)
    else:
        trajectory = OEField("GMXTrajectory_OMD", Types.String, meta=_metaHidden)
        gmx_restart = OEField("GMXRestart_OMD", Types.String, meta=_metaHidden)


class Gromacs:
    gmx_tpr_prefix = "gmx_tpr"

    gmx_mdp_prefix = "gmx_mdp"

    gmx_restart_prefix = "gmx_restart"

    gmx_run_prefix = "gmx_run"

    gmx_traj_dir_prefix = "gmx_traj"

