## Orion MD new Data Record v0.0

---------------------------

Following our meeting during the floe school, We have sketched possible ideas 
on how we would like to implement MD in Orion by using the new OpenEye Data 
Record. We are proposing to use two main OE Records. A first record named 
**MDSystemRecord** will carry info related to the system topology as *OEMol* 
and the synchronized *Parmed Structure* that will hold info related to the 
MD system state (positions and velocities) in addition to the force field:  

|**MDSystemRecord**|  
| ---------------- |  
|*OEMol*    Type: Chem.Mol|  
|*Parmed Structure*  Type: Custom|  

Each MD floe cube will have as input a second data record made of a list of 
**MDStageRecord** and other info if available. The **MDStageRecord** will have 
info related to the *MD Trajectory* in HDF5 format, *Log Data*, 
an **MDSystemRecord** and finally the *Stage Name* (we have to agree on these names)

|**MDStageRecord**|
| --------------- |
|*Stage Name*  Type: String|  
|*MD Trajectory* (HDF5)   Type: Custom|  
|*Log Data*   Type: Text File|  
|**MDSystemRecord**   Type: Custom|  

The list of **MDStageRecord** is an ordered list so we can easily 
retrieve for example the last simulation to carry forward. We did the best we could 
to remember our meeting discussion and feel free to add anything that is relevant 
before I’ll start to implement. Refer to the following figure for an overview.  

[![Data Record MD](https://github.com/oess/openmm_orion/tree/gcalabro_data_record/MD_Data_Record/images/Plan_MD_DataRecord.png)](https://github.com/oess/openmm_orion/tree/gcalabro_data_record/MD_Data_Record/images/Plan_MD_DataRecord.png)

