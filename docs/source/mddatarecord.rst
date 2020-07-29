.. |A|         replace:: Å

#############
MD DataRecord API
#############

These are tutorials for how to use specific floes.

MD DataRecord API overview
==========================

Molecular Dynamics (MD) simulations are notoriously time consuming and
computational resources demanding. In terms of data consuming many information
need to be stored and retrieved such as atomic coordinates, atomic velocities,
forces and energies to cite only few. Extracting and managing these data
is crucial. The MD DataRecord API is an attempt to simplify the storing and
accessing to it. In the OpenEye Datarecord API data is exchanged between
OpenEye cubes in form of records where POD data, custom objects, json data etc.
can be stored and retrieved by using the associated field names and types.
The MD Datarecord API is built on the top of the OpenEye Datarecord API
formatting the record content data produced during the MD runs and providing
getters and setters functions to access to it.

MD Data Record view
-------------------
The MD data produced along MD runs is structured as follow in a md record:

    * the md record contains a sub-record named “md-stages” where md-stage information is saved;

    * each md-stage is a record itself with a name, type, log data, system topology and MD State info.
      The latter is a custom object that stores data such, atomic positions, velocities and box information
      related to the last md frame recorded along a MD run;

    * the md record at the top level can contain a Parmed object used to carry on the system parametrization
    * the md record at the top level can contain a MDComponents object used to carry info related to the different
      system parts such as ligand, protein, cofactors etc.;
    * finally, other info can be present at the top level such as the starting ligand and protein with their names,
      a unique identifier such as the flask id and the flask name and also minor other information

The following diagram shows the structure of the md record.

.. figure_MDRecord:

.. figure:: ./images/MDDataRecord.png
   :width: 600px
   :align: center
   :alt: Structure of the MD Record

   **Structure of the MD Record**

In order to use the MD Datarecord API, the MDOrion package must be installed. The API has been designed
to transparently work locally and in Orion.

Code snippets
-------------------

The following code snippets give an idea on how to use the API.

.. warning::

    In the following examples **record** is an *OpenEye Datarecord* produced
    running the MD floes

.. code:: python

    from MDOrion.Standards.mdrecord import MDDataRecord

    # To use the MD Datarecord API an OERecord is converted into a MD DataRecord
    md_record = MDDataRecord(record)

    # At this point getters and setters can be used to
    # extract info from the record

    # MD Stage name available
    stage_names = md_record.get_stages_names

    # Get a MD Record Stage
    md_stage_production = md_record.get_stage_by_name(stage_names[2])

    # Extract the logging info from a stage
    info_stage = md_record.get_stage_info(stg_name=stage_names[2])

    # Extract the MD State from a stage
    md_set_up_state = md_record.get_stage_state(stg_name=stage_names[0])

    # Extract the Parmed Structure from the md record and synchronize the
    # positions, velocities and box data to the stage state
    pmd_structure = md_record.get_parmed(sync_stage_name=stage_names[2])

    # Extract the OEMol system flask from the record and its title
    flask = md_record.get_flask
    flask_title = md_record.get_title

    # Extract the trajectory file name associated with a stage. In this
    # case the stage trajectory in unpacked and the trajectory name
    # can be used in any md analysis pkg to be loaded
    trj_name = md_record.get_stage_trajectory(stg_name=stage_names[2])

    # Add a new Stage to the md_stages record
    md_record.add_new_stage(stage_name="New_Stage",
                            stage_type="NPT",
                            topology=flask,
                            mdstate=md_set_up_state,
                            data_fn="test.tar.gz")
